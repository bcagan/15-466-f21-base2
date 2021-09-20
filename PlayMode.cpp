#include "PlayMode.hpp"

#include "LitColorTextureProgram.hpp"

#include "DrawLines.hpp"
#include "Mesh.hpp"
#include "Load.hpp"
#include "gl_errors.hpp"
#include "data_path.hpp"

#include <glm/gtc/type_ptr.hpp>

#include<chrono>
#include <random>
#include <iostream>
#include <assert.h>

#define ROT_FACTOR 0.7f
#define PI  3.14159265f
#define ROT_LIMIT 45.0f / 360.0f * 2.f * PI
#define DEGREE_CONVERT 1.0f / 360.0f * 2.f * PI


std::array<std::array<float,4>, 3> bgCols;


GLuint tree_meshes_for_lit_color_texture_program = 0;
Load< MeshBuffer >tree_meshes(LoadTagDefault, []() -> MeshBuffer const * {
	MeshBuffer const *ret = new MeshBuffer(data_path("tree-world.pnct"));
	tree_meshes_for_lit_color_texture_program = ret->make_vao_for_program(lit_color_texture_program->program);
	return ret;
});

Load< Scene > tree_scene(LoadTagDefault, []() -> Scene const * {
	return new Scene(data_path("tree-world.scene"), [&](Scene& scene, Scene::Transform* transform, std::string const& mesh_name) {
		Mesh const& mesh = tree_meshes->lookup(mesh_name);

		scene.drawables.emplace_back(transform);
		Scene::Drawable& drawable = scene.drawables.back();

		drawable.pipeline = lit_color_texture_program_pipeline;

		drawable.pipeline.vao = tree_meshes_for_lit_color_texture_program;
		drawable.pipeline.type = mesh.type;
		drawable.pipeline.start = mesh.start;
		drawable.pipeline.count = mesh.count;
		drawable.bbox.min = mesh.min;
		drawable.bbox.max = mesh.max;

	});
});

//To Do: //Camera is broken if angled initially at all... why?
//Canopy rotate near upper branch rotate limits

PlayMode::PlayMode() : scene(*tree_scene) {

	//Setting up game state variables
	playing = 1;
	level = 1;
	points = 0;
	speedFactor = 1.0f;

	//Creating bg  colors
	float col1[4] = { 0.6f, 0.8f, 1.0f, 1.0f, };
	float col2[4] = { 0.9f, 0.5f, 0.3f, 1.0f, };
	float col3[4] = { 0.7f, 0.25f, 0.15f, 1.0f, };
	bgCols[0][0] = col1[0]; bgCols[0][1] = col1[1]; bgCols[0][2] = col1[2]; bgCols[0][3] = col1[3];
	bgCols[1][0] = col2[0]; bgCols[1][1] = col2[1]; bgCols[1][2] = col2[2]; bgCols[1][3] = col2[3];
	bgCols[2][0] = col3[0]; bgCols[2][1] = col3[1]; bgCols[2][2] = col3[2]; bgCols[2][3] = col3[3];

	//Set up list of paraemters for whether apples should be drawn on tree
	for (size_t appleInd = 0; appleInd < 14; appleInd++) {
		scene.appleBools[appleInd] = false;
	}
	//get pointers to transforms for convenience:
	for (auto& transform : scene.transforms) {
		if (transform.name == "Branch") branch = &transform;
		else if (transform.name == "Canopy") { //Canopy needs to store the height in order to detect closest fruit
			canopy = &transform;
			canopyHeight = transform.position.z;
			canopyDist = glm::vec3(0.0f, 0.0f, canopyHeight);
		}
		else if (transform.name == "Island") island = &transform;
		else if (transform.name == "Fruit0") {  //All fruit need to also store their world positions for later calculations
			fruit0 = &transform;
			fruitPos[0] = transform.make_local_to_world() * glm::vec4(0.0f, 0.0f, 0.0, 1.0f);
		}
		else if (transform.name == "Fruit1") {
			fruit1 = &transform;
			fruitPos[1] = transform.make_local_to_world() * glm::vec4(0.0f, 0.0f, 0.0, 1.0f);
		}
		else if (transform.name == "Fruit2") { 
			fruit2 = &transform;
			fruitPos[2] = transform.make_local_to_world()* glm::vec4(0.0f, 0.0f,0.0, 1.0f);
		}
		else if (transform.name == "Fruit3") {
			fruit3 = &transform;
			fruitPos[3] = transform.make_local_to_world() * glm::vec4(0.0f, 0.0f, 0.0, 1.0f);
		}
		else if (transform.name == "Sphere") {
			sphere = &transform;
		}
	}
	if (branch == nullptr) throw std::runtime_error("Branch not found.");
	if (canopy == nullptr) throw std::runtime_error("Canopy not found.");
	if (island == nullptr) throw std::runtime_error("Island not found.");
	if (fruit0 == nullptr) throw std::runtime_error("Fruit0 not found.");
	if (fruit1 == nullptr) throw std::runtime_error("Fruit1 not found.");
	if (fruit2 == nullptr) throw std::runtime_error("Fruit2 not found.");
	if (fruit3 == nullptr) throw std::runtime_error("Fruit3 not found.");
	if (sphere == nullptr) throw std::runtime_error("Sphere not found.");
	//One of these is not found

	branch_rotation = branch->rotation;
	canopy_rotation = canopy->rotation;
	island_rotation = island->rotation;
	fruit0_rotation = fruit0->rotation;
	fruit1_rotation = fruit1->rotation;
	fruit2_rotation = fruit2->rotation;
	fruit3_rotation = fruit3->rotation;
	sphere_rotation = sphere->rotation;


	//get pointer to camera for convenience:
	if (scene.cameras.size() != 1) throw std::runtime_error("Expecting scene to have exactly one camera, but it has " + std::to_string(scene.cameras.size()));
	camera = &scene.cameras.front();
}

bool PlayMode::isPlaying() {
	return playing != 0;
}

PlayMode::~PlayMode() {
}

void PlayMode::levelUp() { 
	if (points >= 10){ //If level up, despawn all apples and iterate level
		for (size_t whichApple = 0; whichApple < 14; whichApple++) {
			scene.appleBools[whichApple] = false;
		}
		timer = 0.f; //Restart level state
		points = 0;
		level++; //Go to next level
		if (level == 4) {
			level = 1;
			speedFactor *= 1.5; //Ever 3 levels the game resets but at a higher speed
			endTime = 70.0;
		}
		else {
			endTime = 60.0;
		}
	}
}

bool PlayMode::handle_event(SDL_Event const &evt, glm::uvec2 const &window_size) {

	auto vec3Str = [this](glm::vec3 vec) {
		std::string retStr = std::string("(").append(std::to_string(vec.x));
		retStr.append(std::string(","));
		retStr.append(std::to_string(vec.y));
		retStr.append(std::string(","));
		retStr.append(std::to_string(vec.z));
		retStr.append(std::string(")"));
		return retStr;
	};

	auto quatStr = [this](glm::quat quat) {
		std::string retStr = std::string("(").append(std::to_string(quat.w));
		retStr.append(std::string(","));
		retStr.append(std::to_string(quat.x));
		retStr.append(std::string(","));
		retStr.append(std::to_string(quat.y));
		retStr.append(std::string(","));
		retStr.append(std::to_string(quat.z));
		retStr.append(std::string(")"));
		return retStr;
	};

	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.sym == SDLK_ESCAPE) {
			SDL_SetRelativeMouseMode(SDL_FALSE);
			return true;
		} else if (evt.key.keysym.sym == SDLK_a) {
			left.downs += 1;
			left.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.downs += 1;
			right.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.downs += 1;
			up.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.downs += 1;
			down.pressed = true;
			return true;
		}
	} else if (evt.type == SDL_KEYUP) {
		if (evt.key.keysym.sym == SDLK_a) {
			left.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.pressed = false;
			return true;
		}
	} else if (evt.type == SDL_MOUSEBUTTONDOWN) {
		if (SDL_GetRelativeMouseMode() == SDL_FALSE) {
			SDL_SetRelativeMouseMode(SDL_TRUE);
			return true;
		}
	} else if (evt.type == SDL_MOUSEMOTION) {
		if (SDL_GetRelativeMouseMode() == SDL_TRUE) {
			glm::vec2 motion = glm::vec2(
				evt.motion.xrel / float(window_size.y),
				-evt.motion.yrel / float(window_size.y)
			);
			camera->transform->rotation = glm::normalize(
				camera->transform->rotation
				* glm::angleAxis(-motion.x * camera->fovy, glm::vec3(0.0f, 1.0f, 0.0f)));
			glm::vec3 noZPos = camera->transform->position;/*
			cameraRotAngle -= motion.x * camera->fovy;
			camera->transform->rotation = glm::normalize(glm::angleAxis(cameraRotAngle, glm::vec3(0.0f, 0.0f, 1.0f)));*/
			glm::vec3 zPos = camera->transform->position * glm::vec3(0.0,0.0,1.0);
			noZPos.z = 0.0;
			camera->transform->position = 
				glm::mat3_cast(glm::angleAxis(-motion.x*camera->fovy, glm::vec3(0.0f, 0.0f, 1.0f))) * noZPos + zPos;
			return true;
		}
	}

	return false;
}

void PlayMode::update(float elapsed) {
	
	auto rotateIsland = [this]() { //Animates to rotate the island in stage 3
		if ((float)(((int)timer + 10) % 20) < 5) { //Every 20 seconds, animate for 5 seconds, starting at timer = 10
			if (timer - (float)(((int)timer + 10) / 20 - 1)*20.f <= 10.02f) { //Coin flip to decide axis
				//Change so not new randome engien each time
				unsigned seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count(); //Creating seed,
				std::uniform_int_distribution<uint8_t> flip(0, 1);
				whichAxis = (bool) flip(std::default_random_engine(seed));
			}
			float delta = timer - (float)(((int)timer + 10) / 20 -1) * 20.f - 10.f; //Offsets to get delta to 0 instead of 10 for angle calcs
			if (delta >= 2.5f) delta -= 2.5f;
			//Set to loop twice

			//Animation:
			delta /= 2.5f; 
			delta *= (2.f * PI);
			float theta = sin(delta) * 10.0f * DEGREE_CONVERT; 
			if (whichAxis) islandRotAngles = glm::vec2(theta, 0.0f);
			else islandRotAngles = glm::vec2(0.0f, theta); //Smoothly transition between angles of [-10.0,10.0] degrees

		}
		else { //Default resets angles to be safe
			islandRotAngles = glm::vec2(0.0, 0.0);
		}


		//Sets island rotation
		island_rotation = glm::normalize(glm::angleAxis(islandRotAngles.x, glm::vec3(0.0f, 1.0f, 0.0f)));
		island_rotation = glm::normalize(glm::angleAxis(islandRotAngles.y, glm::vec3(1.0f, 0.0f, 0.0f)) * island_rotation);
		island->rotation = island_rotation;

	};

	auto updateFruit = [this](float delta, float timeToDespawn) { //Animates fruit when they are ready to despawn
		if (timeToDespawn >= delta && timeToDespawn - delta <= 0.66667) {
			float theta = 3.f/2.f*16.f * PI * (timeToDespawn - delta); //Updates all fruit heights for 2/3 of a second
			//to rapidly move up and down (excitedly, one could say)
			if (theta >= 8.f * PI) theta *= -2.f;
			float fruitHeightOld = fruitHeightOff;
			fruitHeightOff = (cos(theta + PI) + 1.f) / 2.f; //Speeds up halfway through
			fruit0->position.z += (fruitHeightOff - fruitHeightOld);
			fruit1->position.z += (fruitHeightOff - fruitHeightOld);
			fruit2->position.z += (fruitHeightOff - fruitHeightOld);
			fruit3->position.z += (fruitHeightOff - fruitHeightOld);
		}
	};

	auto spawnApples = [this]() { 
		//Change so not new randome engien each time
		unsigned seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count(); //Creating seed,
		std::uniform_int_distribution<uint8_t> dist(10, 13);
		uint8_t firstPick = dist(std::default_random_engine(seed));
		uint8_t secondPick = 11;
		std::uniform_int_distribution<uint8_t> dist2(11, 13); //Randomly generate which apples is the first second (from 4 and 3)
		if(level == 1){
			switch (firstPick) { //Each case either sets the second apple, or does so and handles collisions with first pick
			case(10):
				secondPick = dist2(std::default_random_engine(seed));
				break;
			case(11):
				secondPick = dist2(std::default_random_engine(seed));
				if (secondPick == 11) secondPick = 10;
				break;
			case(12):
				secondPick = dist2(std::default_random_engine(seed));
				if (secondPick == 12) secondPick = 10;
				break;
			default:
				uint8_t secondPick = dist2(std::default_random_engine(seed)) - 1;
			}
			scene.appleBools[secondPick] = true;
		}
		scene.appleBools[firstPick] = true;
	};

	auto despawnApples = [this]() { //Sets the fruit to all not be spawned
		scene.appleBools[10] = false;
		scene.appleBools[11] = false;
		scene.appleBools[12] = false;
		scene.appleBools[13] = false;
	};

	auto detectPoint = [this]() { //Detects if an apple is spawned, if the tree is rotated close to a block,
		// and if said block contains one of the apples that is spawned

		auto updateApples = [this]() { //Randomly selects the next open apple on the canopy to render

			//Change so not new randome engien each time
			unsigned seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count(); //Creating seed,
			std::uniform_int_distribution<uint8_t> dist(0, 9);
			uint8_t firstPick = dist(std::default_random_engine(seed));
			uint8_t truePick = firstPick;

			do {
				if (scene.appleBools[truePick]) {
					truePick++;
					if (truePick >= 10) truePick = 0;
				}
				else {
					scene.appleBools[truePick] = true;
					return;
				}
			} while (firstPick != truePick);
			throw std::runtime_error("Error spawning picked-up apple\n"); //All apples are already displayed
			//(Should never occur)
		};

		auto vecNorm = [this](glm::vec3 thisVec) {
			return (float)sqrt(thisVec.x * thisVec.x + thisVec.y * thisVec.y + thisVec.z * thisVec.z);
		};

		if (abs(branchRotAngles.x + islandRotAngles.x) >= 0.95f * ROT_LIMIT ||
			abs(branchRotAngles.y + islandRotAngles.y) >= 0.95f * ROT_LIMIT) { //If nearly fully rotated
			std::array<float, 4> appleDists;
			for (size_t whichApple = 0; whichApple < 4; whichApple++) { //For each apple get distance to canopy center
				float dist = vecNorm(branch->make_local_to_world()*glm::vec4(canopyDist,1.0) - fruitPos[whichApple]);
				appleDists[whichApple] = dist;
			}
			float minDist = appleDists[0];
			size_t minInd = 0; //Store the index and distance of closest fruit
			for (size_t whichApple = 1; whichApple < 4; whichApple++) {
				if (appleDists[whichApple] < minDist) {
					minDist = appleDists[whichApple];
					minInd = whichApple;
				}
			}
			if (scene.appleBools[minInd + 10] && abs(branchRotAngles.x +islandRotAngles.x - branchRotAngles.y - +islandRotAngles.y) >= ROT_FACTOR/2) {
				scene.appleBools[minInd + 10] = false; //Only count as a new point if the closest apple has spawned
				points++;
				updateApples();
			}
		}
	};

	//Update spawn/despawn timer
	if (spawnTimer > 0) spawnTimer += speedFactor * elapsed;
	if (despawnTimer > 0) despawnTimer += speedFactor * elapsed;
	//Check to see if new apples should be generated
	if (spawnTimer >= timeToSpawn) {
		spawnApples();
		despawnTimer = spawnTimer - timeToSpawn + 0.0001f; //Error added to allow 0 to be an "off" state
		spawnTimer = 0; //Reset timer and move remainder to other timer
	}
	//Check to see if apples should be despawned
	if (despawnTimer >= timeToDespawn) {
		despawnApples();
		spawnTimer = despawnTimer - timeToDespawn + 0.0001f;
		despawnTimer = 0;
	}
	//If apples are present, check to see if points should be earned
	if (despawnTimer >= 0 && despawnTimer < timeToDespawn) {
		detectPoint();
	}
	//Animation for when despawn is approaching
	updateFruit(despawnTimer, timeToDespawn);

	//Checks game state
	timer += speedFactor * elapsed;
	if (timer >= endTime) { //game over
		playing = 0;
		std::cout << "game over\n";
	}
	levelUp();
	if (level == 3) rotateIsland(); //Level 3 adds island rotations every 20 seconds

	//move camera:
	{

		//combine inputs into a move:
		constexpr float PlayerSpeed = 30.0f;
		glm::vec2 move = glm::vec2(0.0f);
		if (left.pressed && !right.pressed) branchRotAngles.x = branchRotAngles.x - speedFactor * elapsed*ROT_FACTOR;
		if (!left.pressed && right.pressed)  branchRotAngles.x = branchRotAngles.x + speedFactor * elapsed * ROT_FACTOR;
		if (down.pressed && !up.pressed)  branchRotAngles.y = branchRotAngles.y + speedFactor * elapsed * ROT_FACTOR;
		if (!down.pressed && up.pressed) branchRotAngles.y = branchRotAngles.y - speedFactor * elapsed * ROT_FACTOR;
		if (branchRotAngles.x + islandRotAngles.x >= ROT_LIMIT) {
			branchRotAngles.x = ROT_LIMIT - islandRotAngles.x;
			assert(branchRotAngles.x + islandRotAngles.x <= ROT_LIMIT + 0.0001f);
		}
		else if (branchRotAngles.x + islandRotAngles.x <= -ROT_LIMIT) branchRotAngles.x = -ROT_LIMIT - islandRotAngles.x;
		if (branchRotAngles.y + islandRotAngles.y >= ROT_LIMIT) branchRotAngles.y = ROT_LIMIT - islandRotAngles.y;
		else if (branchRotAngles.y + islandRotAngles.y <= -ROT_LIMIT) branchRotAngles.y = -ROT_LIMIT - islandRotAngles.y;

		//Update proper branch rotation
		branch_rotation = glm::normalize(glm::angleAxis(branchRotAngles.x, glm::vec3(0.0f, 1.0f, 0.0f)));
		branch_rotation = glm::normalize(glm::angleAxis(branchRotAngles.y, glm::vec3(1.0f, 0.0f, 0.0f))* branch_rotation);
		branch->rotation = branch_rotation;
		sphere_rotation = branch_rotation;
		sphere->rotation = sphere_rotation;

		//make it so that moving diagonally doesn't go faster:
		if (move != glm::vec2(0.0f)) move = glm::normalize(move) * PlayerSpeed * elapsed;

		glm::mat4x3 frame = camera->transform->make_local_to_parent();
		glm::vec3 right = frame[0];
		glm::vec3 forward = -frame[2];

		camera->transform->position += move.x * right + move.y * forward;
	}

	//reset button press counters:
	left.downs = 0;
	right.downs = 0;
	up.downs = 0;
	down.downs = 0;
}

void PlayMode::draw(glm::uvec2 const &drawable_size) {
	//update camera aspect ratio for drawable:
	camera->aspect = float(drawable_size.x) / float(drawable_size.y);

	glUseProgram(lit_color_texture_program->program);
	glUniform1i(lit_color_texture_program->LIGHT_TYPE_int, 1);
	glUniform3fv(lit_color_texture_program->LIGHT_DIRECTION_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 0.0f,-1.0f)));
	glUniform3fv(lit_color_texture_program->LIGHT_ENERGY_vec3, 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 0.95f)));
	glUseProgram(0);

	std::array<float,4> curCol = bgCols[level - 1];
	glClearColor(curCol[0],curCol[1],curCol[2],curCol[3]);
	glClearDepth(1.0f); //1.0 is actually the default value to clear the depth buffer to, but FYI you can change it.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS); //this is the default depth comparison function, but FYI you can change it.

	GL_ERRORS(); //print any errors produced by this setup code

	scene.draw(*camera);

	{ //use DrawLines to overlay some text:
		glDisable(GL_DEPTH_TEST);
		float aspect = float(drawable_size.x) / float(drawable_size.y);
		DrawLines lines(glm::mat4(
			1.0f / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		));

		constexpr float H = 0.09f;
		std::string timerStr = (std::string("; Time left ")).append(std::to_string(int(endTime - timer)));
		std::string useStr = (std::string("Mouse motion rotates camera; WASD moves; escape ungrabs mouse")).append(timerStr);
		lines.draw_text(useStr.c_str(),
			glm::vec3(-aspect + 0.1f * H, -1.0 + 0.1f * H, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0x00, 0x00, 0x00, 0x00));
		float ofs = 2.0f / drawable_size.y;
		lines.draw_text(useStr.c_str(),
			glm::vec3(-aspect + 0.1f * H + ofs, -1.0 + + 0.1f * H + ofs, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0xff, 0xff, 0xff, 0x00));
	}
}
