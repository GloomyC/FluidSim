
#include <iostream>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <random>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <math.h>
#include "SimParts.h"
#include <thread>

using namespace sf;
using namespace std;

Color pressure(float p) {
	float k = 100.0f;
	if (p * k <= -255) return Color::Black;
	if (p * k <= 0) return Color::Color(0, 0, 255 + (int(p * k)), 255);;
	if (p * k <= 255) return Color::Color(int(p * k), 0, 255 - (int(p * k)), 255);
	return Color::Red;
}

template <class T>
Color randColor(T i) {
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(0, 255);

	return Color(dis(gen), dis(gen), dis(gen), 255);
};

void main()
{
	Eigen::initParallel();
	cout << "asdf" << thread::hardware_concurrency << endl;
	int sizeX = 10;
	int sizeY = 10;

	RenderWindow window{ VideoMode{800,800}, "Tutorial" };
	window.setFramerateLimit(1);
	Event event;
	CellBoard<float> board(Vector2i(800, 800), Vector2i(sizeX, sizeY));
	board.setColorF(pressure);

	FluidSim sim(10, 1, sizeX, sizeY, 1.0f, 1.0f);
	board.setBoard(sim.getRaw());

	Vector2i prevclick;
	bool flag = true;
	Matrix<float> vels(sizeX, sizeY, 0.0f, 0.0f);
	while (window.isOpen()) {

		window.clear(Color::Black);
		if (window.pollEvent(event)) {
			if (event.type == Event::Closed) {
				window.close();
				break;
			}
			if (Mouse::isButtonPressed(Mouse::Button::Left)) {
				Vector2i cellPos = board.getCellCoordAtPos(Mouse::getPosition(window));
				sim.clearvX(prevclick.x, prevclick.y);
				prevclick = cellPos;
				sim.wiatrak(prevclick.x, prevclick.y);
			}
			if (event.type == Event::MouseButtonPressed && event.mouseButton.button == Mouse::Right) {
				Vector2i cellPos = board.getCellCoordAtPos(Mouse::getPosition(window));
				flag = !flag;

			}


		}
		sim.sim();

		for (int i = 0; i < sizeX; i++) {
			for (int j = 0; j < sizeY; j++) {
				vels(i, j) = sim.interpVx(i, j);
			}
		}
		if (flag) {
			board.setBoard(sim.getRaw());
		}
		else {
			board.setBoard(vels.getPtr());
		}
		board.refreshLook();
		window.draw(board);
		window.display();
	}

}

