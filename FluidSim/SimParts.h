
#ifndef SIM_PARTS
#define SIM_PARTS

#include <iostream>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <random>
#include <Eigen/Sparse>
#include <math.h>
#include <thread>

using namespace std;
using namespace sf;

static const int thread_num = 4;

template <class T>
class CellBoard : public Drawable, public Transformable {
private:
	T* values;
	uint8_t* pixels;
	Texture texture;
	Color(*valtoColor)(T);
	VertexArray shape;
	Vector2i cellsSize;
	Vector2i size;

	uint8_t* getPixel(int x, int y) {
		return &pixels[(y * cellsSize.x + x) * 4];
	}



public:
	//constructors
	CellBoard(Vector2i size, Vector2i cells) {
		cellsSize = cells;
		this->size = size;
		pixels = new uint8_t[(cells.x * cells.y) * 4];
		texture.create(cells.x, cells.y);

		shape = VertexArray(Quads, 4);
		shape[0].position = Vector2f(0, 0);
		shape[1].position = Vector2f(size.x, 0);
		shape[2].position = Vector2f(size.x, size.y);
		shape[3].position = Vector2f(0, size.y);
		shape[0].texCoords = Vector2f(0, 0);
		shape[1].texCoords = Vector2f(cells.x, 0);
		shape[2].texCoords = Vector2f(cells.x, cells.y);
		shape[3].texCoords = Vector2f(0, cells.y);


		valtoColor = [](T i) {return Color::Black; };

	}
	~CellBoard() {
		delete[] pixels;
	}

	//getters
	T* operator()(int x, int y) {
		if (isInBounds(x, y)) return &values[y * cellsSize.x + x];
		return NULL;
	}

	Vector2i getCellsSize() {
		return cellsSize;
	}

	Vector2i getSize() {
		return size;
	}

	T* getBoard() {
		return (*this)(0, 0);
	}

	Vector2i getCellCoordAtPos(Vector2i pos) {
		return Vector2i(pos.x / (size.x / cellsSize.x), pos.y / (size.y / cellsSize.y));
	}

	//modify

	void setBoard(T* newB) {
		this->values = newB;
	}

	void fillAll(T val) {
		for (int i = 0; i < cellsSize.x * cellsSize.y; i++) {
			values[i] = val;
		}
	}

	void setColorF(Color(*function)(T)) {
		valtoColor = function;
	}


	bool isInBounds(int x, int y) {
		if (x < 0 || x >= cellsSize.x || y < 0 || y >= cellsSize.y) return false;
		return true;
	}

	void refreshLook() {
		for (int y = 0; y < cellsSize.y; y++) {
			for (int x = 0; x < cellsSize.x; x++) {
				Color newC = valtoColor(*(*this)(x, y));
				uint8_t* pixel = getPixel(x, y);

				pixel[0] = newC.r;
				pixel[1] = newC.g;
				pixel[2] = newC.b;
				pixel[3] = newC.a;
			}
		}
		texture.update(pixels);
	}

	virtual void draw(RenderTarget& target, RenderStates states) const {

		states.transform *= getTransform();
		states.texture = &texture;
		target.draw(shape, states);

	}





};

template<class T>
class Matrix {
private:
	int x;
	int y;
	T* vals;
	T outOfBoundsVal;
public:
	Matrix() {
	};
	Matrix(int x, int y, T init, T outOfBounds) {
		this->x = x;
		this->y = y;
		this->outOfBoundsVal = outOfBounds;
		vals = new T[x * y];
		for (int i = 0; i < x * y; i++) {
			vals[i] = init;
		}
	}
	~Matrix() {
		delete[] vals;
	}
	Matrix(const Matrix& m) {

		x = m.x;
		y = m.y;
		vals = new T[x * y];
		for (int i = 0; i < x * y; i++) {
			vals[i] = m.vals[i];
		}
		outOfBoundsVal = m.outOfBoundsVal;
	}
	Matrix& operator=(const Matrix& m) {
		if (this == &m) {
			return *this;
		}

		this->x = m.x;
		this->y = m.y;
		this->vals = new T[x * y];
		for (int i = 0; i < x * y; i++) {
			vals[i] = m.vals[i];
		}
		this->outOfBoundsVal = m.outOfBoundsVal;
		return *this;
	}

	int X() { return x; }
	int Y() { return y; }

	T& operator()(int x, int y) {
		if (inBounds(x, y)) return vals[y * this->x + x];
		return outOfBoundsVal;
	}
	bool inBounds(int x, int y) {
		if (x < 0 || x >= this->x || y < 0 || y >= this->y) return false;
		return true;
	}
	T* getPtr() { return vals; }

	void setPtr(float* ptr, bool destroy) {
		if (destroy) delete[] vals;
		vals = ptr;
	}

};

class FluidSim {
public:
	int sizeX, sizeY;
	float dx;
	float dt;
	float density;
	float visc;
	int maxIter;

	Matrix<float> vx0;
	Matrix<float> vy0;
	Matrix<float> p0;

	Matrix<float> vx;
	Matrix<float> vy;
	Matrix<float> p;


	FluidSim(float density, float visc, int sizeX, int sizeY, float dt = 0.1f, float dx = 1.0f) {
		this->dt = dt;
		this->dx = dx;
		this->density = density;
		this->visc = visc;
		this->sizeX = sizeX;
		this->sizeY = sizeY;
		maxIter = 10;

		//0 are current state, 
		//without nmr are next state
		vx = Matrix<float>(sizeX + 1, sizeY, 0.0f, 0.0f);
		vx0 = Matrix<float>(sizeX + 1, sizeY, 0.0f, 0.0f);
		vy = Matrix<float>(sizeX, sizeY + 1, 0.0f, 0.0f);
		vy0 = Matrix<float>(sizeX, sizeY + 1, 0.0f, 0.0f);
		p0 = Matrix<float>(sizeX, sizeY, 0.0f, 0.0f);
	}
	~FluidSim() {
	}

	float divergence(int x, int y) {
		return ((vx0.inBounds(x, y) ? vx0(x, y) : 0.0f)
			- (vx0.inBounds(x + 1, y) ? vx0(x + 1, y) : 0.0f)
			+ (vy0.inBounds(x, y) ? vy0(x, y) : 0.0f)
			- (vy0.inBounds(x, y + 1) ? vy0(x, y + 1) : 0.0f)) / dx;
	}

	//calculate pressure so that fluid stays incompressible
	void project() {
		//prepare b vector
		Eigen::VectorXf b(sizeX * sizeY);

		thread threads[thread_num];
		int launched_threads = 0;
		for (int j = 0; j < sizeY; ++j) {

			//lambda setting divergences in given row
			auto diverRow = [](int j, Eigen::VectorXf* target, FluidSim* self)->void {
				for (int i = 0; i < self->sizeX; i++) {
					(*target)[j * self->sizeX + i] = self->divergence(i, j);
					cout << self->divergence(i, j);
				}
			};

			//launch thread if avilable
			if (launched_threads < thread_num) {
				threads[launched_threads] = thread(diverRow, j, &b, this);
				launched_threads++;
			}
			else {
				//wait for all threads if none are avilable
				for (int n = 0; n < thread_num; n++) {
					threads[n].join();
				}
				launched_threads = 0;
			}


		}
		//wait for threads
		for (int n = 0; n < launched_threads; n++) {
			threads[n].join();
		}
		launched_threads = 0;

		cout << "------------------------" << endl;
		cout << b << endl;

		cout << "------------------------" << endl;



		//prepare matrix a
		Eigen::SparseMatrix<float> a(sizeX * sizeY, sizeX * sizeY);
		for (int i = 0; i < sizeX; i++) {
			for (int j = 0; j < sizeY; j++) {
				int count = 0;

				float k = dt * (density * dx * dx);
				if (p0.inBounds(i - 1, j)) {
					count++;
					a.insert(j * sizeX + i, j * sizeX + i - 1) = -k;
				}
				if (p0.inBounds(i + 1, j)) {
					count++;
					a.insert(j * sizeX + i, j * sizeX + i + 1) = -k;
				}
				if (p0.inBounds(i, j - 1)) {
					count++;
					a.insert(j * sizeX + i, (j - 1) * sizeX + i) = -k;
				}
				if (p0.inBounds(i, j + 1)) {
					count++;
					a.insert(j * sizeX + i, (j + 1) * sizeX + i) = -k;
				}
				a.insert(j * sizeX + i, j * sizeX + i) = float(count) * k;
			}
		}

		//solve
		Eigen::ConjugateGradient< Eigen::SparseMatrix<float>> solver;
		solver.compute(a);
		solver.setMaxIterations(maxIter);

		Eigen::VectorXf solution(sizeX * sizeY);
		solution = solver.solve(b);
		cout << solver.iterations() << endl;

		for (int i = 0; i < sizeX; i++) {
			for (int j = 0; j < sizeY; j++) {
				p0(i, j) = solution[j * sizeX + i];
			}
		}


		//update velocities from pressure
		for (int i = 1; i < sizeX; i++) {
			for (int j = 0; j < sizeY; j++) {
				vx(i, j) = vx0(i, j) - dt * (1 / density) * ((p0(i, j) - p0(i - 1, j)) / dx);
			}
		}
		for (int i = 0; i < sizeX; i++) {
			for (int j = 1; j < sizeY; j++) {
				vy(i, j) = vy0(i, j) - dt * (1 / density) * ((p0(i, j) - p0(i, j - 1)) / dx);
			}
		}

		bounds();
		swap();
	}

	//interpolates x part of velocities on board, 
	//center of cell (0,0) is center of coordinates
	float interpVx(float x, float y) {

		//if out of bounds take border value
		if (x < -0.5f) return interpVx(-0.5f, y);
		if (x > float(sizeX) + 0.5f) return interpVx(float(sizeX) + 0.5f, y);
		if (y < 0.0f) return interpVx(x, 0.0f);
		if (y > float(sizeY) - 1.0f) return interpVx(x, float(sizeY) - 1.0f);
		//zero to coordinates of matrix vx0 (grid is staggered)
		x = x + 0.5f;
		//coordinates inside cell
		float epsX = fmod(x, 1.0f);
		float epsY = fmod(y, 1.0f);

		//interpolate from 4 closest data points
		float av0 = (1.0f - epsX) * vx0((int)x, (int)y) + epsX * vx0((int)x + 1, (int)y);
		float av1 = (1.0f - epsX) * vx0((int)x, (int)y + 1) + epsX * vx0((int)x + 1, (int)y + 1);
		return ((1.0f - epsY) * av0) + (epsY * av1);
	}

	//interpolates y part of velocities on board, 
	//center of cell (0,0) is center of coordinates
	float interpVy(float x, float y) {

		//if out of bounds take border value
		if (x < 0.0f) return interpVy(0.0f, y);
		if (x > float(sizeX) - 1.0f) return interpVy(float(sizeX) - 1.0f, y);
		if (y < -0.5f) return interpVy(x, -0.5f);
		if (y > float(sizeY) + 0.5f) return interpVy(x, float(sizeY) + 0.5f);
		//zero to coordinates of matrix vy0 (grid is staggered)
		y = y + 0.5f;
		//coordinates inside cell
		float epsX = fmod(x, 1.0f);
		float epsY = fmod(y, 1.0f);

		//interpolate from 4 closest data points
		float av0 = (1.0f - epsX) * vy0((int)x, (int)y) + epsX * vy0((int)x + 1, (int)y);
		float av1 = (1.0f - epsX) * vy0((int)x, (int)y + 1) + epsX * vy0((int)x + 1, (int)y + 1);
		return ((1.0f - epsY) * av0) + (epsY * av1);
	}

	void advect() {
		//advect vx
		for (int i = 0; i < sizeX + 1; i++) {
			for (int j = 0; j < sizeY; j++) {
				float xMid = (float)i - 0.5f * dt * vx0(i, j);
				float yMid = (float)j - 0.5f * dt * interpVy((float)i - 0.5f, (float)j);

				float xEnd = (float)i - dt * interpVx(xMid, yMid);

				vx(i, j) = xEnd;

			}
		}

		for (int i = 0; i < sizeX; i++) {
			for (int j = 0; j < sizeY + 1; j++) {
				float yMid = (float)j - 0.5f * dt * vy0(i, j);
				float xMid = (float)i - 0.5f * dt * interpVx((float)i, (float)j - 0.5f);

				float yEnd = (float)j - dt * interpVy(xMid, yMid);

				vy(i, j) = yEnd;

			}
		}

		bounds();
		swap();
	}
	void bounds() {
		for (int i = 0; i < sizeX; i++) {
			vy(i, 0) = 0.0f;
			vy(i, sizeY) = 0.0f;
		}
		for (int j = 0; j < sizeY; j++) {
			vy(0, j) = 0.0f;
			vy(sizeX, j) = 0.0f;
		}
	}
	void swap() {
		float* tmp = vx0.getPtr();

		vx0.setPtr(vx.getPtr(), false);
		vx.setPtr(tmp, false);

		tmp = vy0.getPtr();
		vy0.setPtr(vy.getPtr(), false);
		vy.setPtr(tmp, false);
	}

	void sim() {

		//advect();
		project();
	}

	float* getRaw() {
		return p0.getPtr();
	}

	void testVelocity() {
		vx0(5, 5) = 1.0f;
	}

	void testPressure() {
		p0(0, 0) = 200.0f;
	}
	void wiatrak(int x, int y) {
		vx0(x + 1, y) = 100.0f;
	}
	void clearvX(int x, int y) {

		vx0(x + 1, y) = 0.0f;

	}
};

#endif 