#include <random>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/random.hpp>
#include <ctime>
#include <cstdlib>
#include <cmath>

/**
 * The number of particles we want to simulate. Since compute shaders are limited in their
 * work group size, we'll also need to know the work group size so we can find out how many
 * work groups to dispatch
 */
#define NUM_PARTICLES 100

// This MUST match with local_size_x inside cursor.glsl
#define WORK_GROUP_SIZE 10

#define SCREENX 1920
#define SCREENY 1080

GLchar *LoadShader(const std::string &file)
{
	std::ifstream shaderFile;
	long shaderFileLength;

	shaderFile.open(file);

	if (shaderFile.fail())
	{
		throw std::runtime_error("COULD NOT FIND SHADER FILE");
	}

	shaderFile.seekg(0, shaderFile.end);
	shaderFileLength = shaderFile.tellg();
	shaderFile.seekg(0, shaderFile.beg);

	GLchar *shaderCode = new GLchar[shaderFileLength + 1];
	shaderFile.read(shaderCode, shaderFileLength);

	shaderFile.close();

	shaderCode[shaderFileLength] = '\0';

	return shaderCode;
}

int main(int argc, char **argv)
{
	srand(time(0));

	const GLfloat delta_time = 0.01f;

	// Window Setup
	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);

	GLFWwindow *window = glfwCreateWindow(SCREENX, SCREENY, "Particles", nullptr, nullptr);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	gladLoadGL();
	
	// Setup Initial Positions/Velocities
	glm::vec3 positions[NUM_PARTICLES];
	glm::vec3 velocities[NUM_PARTICLES];

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, M_PI);
	std::uniform_real_distribution<double> dis2(0, M_PI);
	std::uniform_real_distribution<double> dis3(1, 10);
	std::uniform_int_distribution<int> dis4(1, 10);

	for (int i = 0; i < NUM_PARTICLES; i++)
	{
		double theta = dis1(gen);
		double phi = dis2(gen);
		double rand_mass = dis3(gen);
		int r = dis4(gen);

		float velx = r * sin(theta) * cos(phi);
		float vely = r * sin(theta) * sin(phi);
		float velz = r * cos(theta);

		positions[i] = glm::vec3(velx, vely, velz);
		velocities[i] = glm::vec3(0.0f, 0.0f, 0.0f);
	}


	GLuint vao, pos, timebuffer, vel;

	glCreateVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glCreateBuffers(1, &pos);
	glCreateBuffers(1, &timebuffer);
	glCreateBuffers(1, &vel);

	glNamedBufferData(vel, sizeof(velocities), velocities, GL_STATIC_DRAW);
	glNamedBufferData(pos, sizeof(positions), positions, GL_STATIC_DRAW);
	glNamedBufferData(timebuffer, sizeof(GLfloat), &delta_time, GL_DYNAMIC_DRAW);

	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, pos);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, vel);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, timebuffer);
	/*
	 * Basic Shader Setup - This is standard shader loading and compilation.
	 */
	const GLchar *vertCode = LoadShader("../shader_files/shader.vert");
	const GLchar *fragCode = LoadShader("../shader_files/shader.frag");
	const GLchar *computeCode = LoadShader("../shader_files/NBodies.glsl");

	GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
	GLuint computeShader = glCreateShader(GL_COMPUTE_SHADER);

	glShaderSource(vertShader, 1, &vertCode, nullptr);
	glShaderSource(fragShader, 1, &fragCode, nullptr);
	glShaderSource(computeShader, 1, &computeCode, nullptr);

	GLchar infolog[512];

	glCompileShader(vertShader);
	glGetShaderInfoLog(vertShader, 512, nullptr, infolog);
	std::cout << infolog << std::endl;

	glCompileShader(fragShader);
	glGetShaderInfoLog(fragShader, 512, nullptr, infolog);
	std::cout << infolog << std::endl;

	glCompileShader(computeShader);
	glGetShaderInfoLog(computeShader, 512, nullptr, infolog);
	std::cout << infolog << std::endl;

	GLuint shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, vertShader);
	glAttachShader(shaderProgram, fragShader);
	glLinkProgram(shaderProgram);
	glDeleteShader(vertShader);
	glDeleteShader(fragShader);

	glGetProgramInfoLog(shaderProgram, 512, nullptr, infolog);

	std::cout << infolog << std::endl;

	GLuint computeProgram = glCreateProgram();
	glAttachShader(computeProgram, computeShader);
	glLinkProgram(computeProgram);
	glDeleteShader(computeShader);

	glGetProgramInfoLog(computeProgram, 512, nullptr, infolog);
	std::cout << infolog << std::endl;

	/**
	 * Setup some basic properties for screen size, background color and point size
	 */
	glViewport(0, 0, SCREENX, SCREENY);
	glClearColor(0.05f, 0.05, 0.05f, 1.0f);
	glEnable(GL_DEPTH_TEST);
	glPointSize(2.0f);

	// Draw Loop
	while (glfwWindowShouldClose(window) == 0)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glfwPollEvents();

		glUseProgram(computeProgram);
		glDispatchCompute(NUM_PARTICLES / WORK_GROUP_SIZE, 1, 1);

		/**
		 * Draw the particles
		 */
		glUseProgram(shaderProgram);
		glDrawArraysInstanced(GL_POINTS, 0, 1, NUM_PARTICLES);

		glfwSwapBuffers(window);
	}

	// OpenGL Shutdown
	glDeleteProgram(shaderProgram);
	glDeleteProgram(computeProgram);
	delete[] vertCode;
	delete[] fragCode;
	delete[] computeCode;

	// Window Shutdown
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}