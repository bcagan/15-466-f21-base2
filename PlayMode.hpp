#include "Mode.hpp"

#include "Scene.hpp"

#include <glm/glm.hpp>

#include <vector>
#include <deque>
#include <array>

struct PlayMode : Mode {
	PlayMode();
	virtual ~PlayMode();

	//functions called by main loop:
	virtual bool handle_event(SDL_Event const &, glm::uvec2 const &window_size) override;
	virtual void update(float elapsed) override;
	virtual void draw(glm::uvec2 const &drawable_size) override;

	//----- game state -----

	//input tracking:
	struct Button {
		uint8_t downs = 0;
		uint8_t pressed = 0;
	} left, right, down, up;

	float points = 0;
	size_t level = 1;

	//Times for spawning apples
	float spawnTimer = 0.000001f;
	float despawnTimer = 0.f;
	float timeToSpawn = 3.f;
	float timeToDespawn = 2.1f;
	float fruitHeightOff = 0.0f;

	//Game timer
	float timer = 0.f;
	float endTime = 60.0f;
	virtual bool isPlaying() override;
	int playing = 1;
	void levelUp();
	float speedFactor = 1.0f;

	//local copy of the game scene (so code can change it during gameplay):
	Scene scene;

	//Camera rot
	float cameraRotAngle = 0.0f;

	//Tree Transforms
	Scene::Transform *branch = nullptr;
	Scene::Transform *canopy = nullptr;
	glm::quat canopy_rotation;
	glm::quat branch_rotation;
	glm::vec2 branchRotAngles = glm::vec2(0.0, 0.0);
	float canopyHeight;

	//Stage Transforms
	//(Fruit will wobble when they first appear, and then start shaking up and down right before they disappear)
	Scene::Transform* fruit0 = nullptr;
	glm::quat fruit0_rotation;
	Scene::Transform* fruit1 = nullptr;
	glm::quat fruit1_rotation;
	Scene::Transform* fruit2 = nullptr;
	glm::quat fruit2_rotation;
	Scene::Transform* fruit3 = nullptr;
	glm::quat fruit3_rotation;
	//Island can rotate itself during stage 3
	Scene::Transform* island = nullptr;
	glm::quat island_rotation;
	std::array<glm::vec3, 4> fruitPos;
	glm::vec3 canopyDist;

	float wobble = 0.0f;
	
	//camera:
	Scene::Camera *camera = nullptr;

};
