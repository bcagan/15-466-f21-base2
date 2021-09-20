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
// 
//Make level 3 (PRIORITY OVER 2) have an island rotate transform (ANIMATION, pick from a few and start at random if none are occuring)
// Include island rotate with LOCK LIMIT ONLY of branch rotate
// 
//Canopy rotate near upper branch rotate limits
//Loop and speed up after level 3 (item spawn, animation speed, and control spped)
//Make asset textures change with each level

PlayMode::PlayMode() : scene(*tree_scene) {

	playing = 1;
	level = 1;
	points = 0;
	float speedFactor = 1.0f;

	auto getBBox = [this](std::string name) {
		for (auto& thisDraw : scene.drawables) {
			if (thisDraw.transform->name.compare(name) == 0) {
				return thisDraw.bbox;
			}
		}
		MeshBBox newBox;
		newBox.min = glm::vec3(0.0f);
		newBox.max = newBox.min;
		std::string errString = std::string("BBox not found not found. Mesh: %s").append(name);
		throw std::runtime_error(errString.c_str());
		return newBox;
	};

	//Set up list of paraemters for whether apples should be drawn on tree
	for (size_t appleInd = 0; appleInd < 14; appleInd++) {
		scene.appleBools[appleInd] = false;
	}
	//get pointers to transforms for convenience:
	for (auto& transform : scene.transforms) {
		if (transform.name == "Branch") branch = &transform;
		else if (transform.name == "Canopy") {
			canopy = &transform;
			MeshBBox canopyBox = getBBox(std::string("Canopy"));
			float distX = abs(canopyBox.max.x - canopyBox.min.x);
			float distY = abs(canopyBox.max.y - canopyBox.min.y);
			canopyHeight = transform.position.z;
			canopyDist = glm::vec3(0.0f, 0.0f, canopyHeight);
		}
		else if (transform.name == "Island") island = &transform;
		else if (transform.name == "Fruit0") { 
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
	}
	if (branch == nullptr) throw std::runtime_error("Branch not found.");
	if (canopy == nullptr) throw std::runtime_error("Canopy not found.");
	if (island == nullptr) throw std::runtime_error("Island not found.");
	if (fruit0 == nullptr) throw std::runtime_error("Fruit0 not found.");
	if (fruit1 == nullptr) throw std::runtime_error("Fruit1 not found.");
	if (fruit2 == nullptr) throw std::runtime_error("Fruit2 not found.");
	if (fruit3 == nullptr) throw std::runtime_error("Fruit3 not found.");
	//One of these is not found

	branch_rotation = branch->rotation;
	canopy_rotation = canopy->rotation;
	island_rotation = island->rotation;
	fruit0_rotation = fruit0->rotation;
	fruit1_rotation = fruit1->rotation;
	fruit2_rotation = fruit2->rotation;
	fruit3_rotation = fruit3->rotation;



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
	for (size_t whichApple = 0; whichApple < 14; whichApple++) {
		scene.appleBools[whichApple] = false;
	}
	level = 2;
	std::cout << "level up!\n";
	timer = 0.f;
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
				glm::mat3_cast(glm::angleAxis(-motion.x * camera->fovy, glm::vec3(0.0f, 0.0f, 1.0f))) * noZPos + zPos;
			return true;
		}
	}

	return false;
}

void PlayMode::update(float elapsed) {
	/*unsigned seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count(); //Creating seed,
	//Found seed function from http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
	std::uniform_real_distribution < double > dist(0.0, 100.0 * (double) maxTop - 100.0*(double)(curGap + minBottom));
	float curTop = (float) dist(std::default_random_engine(seed))/100.f + curGap + minBottom; //Randomly make new gate (top of gate gap)
	dist.reset();
	*/

	auto updateFruit = [this](float delta, float timeToDespawn) {
		if (timeToDespawn >= delta && timeToDespawn - delta <= 0.66667) {
			float theta = 3.f/2.f*16.f * PI * (timeToDespawn - delta);
			if (theta >= 8.f * PI) theta *= -2.f;
			float fruitHeightOld = fruitHeightOff;
			fruitHeightOff = (cos(theta + PI) + 1.f) / 2.f;
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
		std::uniform_int_distribution<uint8_t> dist2(11, 13);
		if(level == 1){
			switch (firstPick) {
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

	auto despawnApples = [this]() {
		scene.appleBools[10] = false;
		scene.appleBools[11] = false;
		scene.appleBools[12] = false;
		scene.appleBools[13] = false;
	};

	auto detectPoint = [this]() {

		auto updateApples = [this]() {

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
			throw std::runtime_error("Error spawning picked-up apple\n");
		};

		auto vecNorm = [this](glm::vec3 thisVec) {
			return (float)sqrt(thisVec.x * thisVec.x + thisVec.y * thisVec.y + thisVec.z * thisVec.z);
		};

		//Only works for one apple
		if (abs(branchRotAngles.x) >= 0.95f * ROT_LIMIT || abs(branchRotAngles.y) >= 0.95f * ROT_LIMIT) { //If nearly fully rotated
			std::array<float, 4> appleDists;
			for (size_t whichApple = 0; whichApple < 4; whichApple++) {
				float dist = vecNorm(branch->make_local_to_world()*glm::vec4(canopyDist,1.0) - fruitPos[whichApple]);
				appleDists[whichApple] = dist;
			}
			float minDist = appleDists[0];
			size_t minInd = 0;
			for (size_t whichApple = 1; whichApple < 4; whichApple++) {
				if (appleDists[whichApple] < minDist) {
					minDist = appleDists[whichApple];
					minInd = whichApple;
				}
			}
			if (scene.appleBools[minInd + 10] && abs(branchRotAngles.x - branchRotAngles.y) >= ROT_FACTOR/2) {
				scene.appleBools[minInd + 10] = false;
				points++;
				updateApples();
			}
		}
	};

	//Update spawn/despawn timer
	if (spawnTimer > 0) spawnTimer += elapsed;
	if (despawnTimer > 0) despawnTimer += elapsed;
	if (spawnTimer >= timeToSpawn) {
		spawnApples();
		despawnTimer = spawnTimer - timeToSpawn + 0.0001f; //Error added to allow 0 to be an "off" state
		spawnTimer = 0; //Reset timer and move remainder to other timer
	}
	if (despawnTimer >= timeToDespawn) {
		despawnApples();
		spawnTimer = despawnTimer - timeToDespawn + 0.0001f;
		despawnTimer = 0;
	}
	if (despawnTimer >= 0 && despawnTimer < timeToDespawn) {
		detectPoint();
	}
	updateFruit(despawnTimer, timeToDespawn);

	timer += elapsed;
	std::cout << "Time left " << (int)(endTime - timer) << std::endl;
	if (timer >= endTime) { //game over
		playing = 0;
		std::cout << "game over\n";
	}
	if (points >= 10 && level == 1) {
		levelUp();
		points = 0;
	}

	//slowly rotates through [0,1):
	/*
	wobble += elapsed / 10.0f;
	wobble -= std::floor(wobble);

	hip->rotation = hip_base_rotation * glm::angleAxis(
		glm::radians(5.0f * std::sin(wobble * 2.0f * float(M_PI))),
		glm::vec3(0.0f, 1.0f, 0.0f)
	);
	upper_leg->rotation = upper_leg_base_rotation * glm::angleAxis(
		glm::radians(7.0f * std::sin(wobble * 2.0f * 2.0f * float(M_PI))),
		glm::vec3(0.0f, 0.0f, 1.0f)
	);
	lower_leg->rotation = lower_leg_base_rotation * glm::angleAxis(
		glm::radians(10.0f * std::sin(wobble * 3.0f * 2.0f * float(M_PI))),
		glm::vec3(0.0f, 0.0f, 1.0f)
	);*/

	//move camera:
	{

		//combine inputs into a move:
		constexpr float PlayerSpeed = 30.0f;
		glm::vec2 move = glm::vec2(0.0f);
		if (left.pressed && !right.pressed) branchRotAngles.x = branchRotAngles.x - elapsed*ROT_FACTOR;
		if (!left.pressed && right.pressed)  branchRotAngles.x = branchRotAngles.x + elapsed * ROT_FACTOR;
		if (down.pressed && !up.pressed)  branchRotAngles.y = branchRotAngles.y + elapsed * ROT_FACTOR;
		if (!down.pressed && up.pressed) branchRotAngles.y = branchRotAngles.y - elapsed * ROT_FACTOR;
		if (branchRotAngles.x >= ROT_LIMIT) branchRotAngles.x = ROT_LIMIT;
		else if (branchRotAngles.x <= -ROT_LIMIT) branchRotAngles.x = -ROT_LIMIT;
		if (branchRotAngles.y >= ROT_LIMIT) branchRotAngles.y = ROT_LIMIT;
		else if (branchRotAngles.y <= -ROT_LIMIT) branchRotAngles.y = -ROT_LIMIT;
		assert(abs(branchRotAngles.x) <= ROT_LIMIT);
		//Update proper branch rotation
		branch_rotation = glm::normalize(glm::angleAxis(branchRotAngles.x, glm::vec3(0.0f, 1.0f, 0.0f)));
		branch_rotation = glm::normalize(glm::angleAxis(branchRotAngles.y, glm::vec3(1.0f, 0.0f, 0.0f))* branch_rotation);
		branch->rotation = branch_rotation;

		//make it so that moving diagonally doesn't go faster:
		if (move != glm::vec2(0.0f)) move = glm::normalize(move) * PlayerSpeed * elapsed;

		glm::mat4x3 frame = camera->transform->make_local_to_parent();
		glm::vec3 right = frame[0];
		//glm::vec3 up = frame[1];
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

	//set up light type and position for lit_color_texture_program:
	// TODO: consider using the Light(s) in the scene to do this
	glUseProgram(lit_color_texture_program->program);
	glUniform1i(lit_color_texture_program->LIGHT_TYPE_int, 1);
	glUniform3fv(lit_color_texture_program->LIGHT_DIRECTION_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 0.0f,-1.0f)));
	glUniform3fv(lit_color_texture_program->LIGHT_ENERGY_vec3, 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 0.95f)));
	glUseProgram(0);

	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
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
		lines.draw_text("Mouse motion rotates camera; WASD moves; escape ungrabs mouse",
			glm::vec3(-aspect + 0.1f * H, -1.0 + 0.1f * H, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0x00, 0x00, 0x00, 0x00));
		float ofs = 2.0f / drawable_size.y;
		lines.draw_text("Mouse motion rotates camera; WASD moves; escape ungrabs mouse",
			glm::vec3(-aspect + 0.1f * H + ofs, -1.0 + + 0.1f * H + ofs, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0xff, 0xff, 0xff, 0x00));
	}
}
