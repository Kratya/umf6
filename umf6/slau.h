#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>
#include <vector>
#include <set>
#include <ctime>

using namespace std;

class SLAU {
public:
	//параметры
	double lyam;
	double omega;
	double hi;
	double sigma;

	int n_x, n_y, n_z; //узлы областей
	int count_x, count_y, count_z; //узлы областей + внутренние

	int max_it = 2000;
	double err = 1e-10;
	int iter;
	double nev;

	vector<double> ax, ay, az; //узлы областей
	vector<double> x, y, z; //узлы-просто все уже узлы
	vector<int> nx, ny, nz;//количество разбиений области
	vector<double> qx, qy, qz;//коэффициент разбиений по области
	vector<int> oxn, oyn, ozn;//еще массивчики соответствий узлов и номера в х,y,z

	int m; //количество глобальных узлов

//	double w = 1;//параметр релаксации

	//краевые
	int ncond1, ncond2;
	vector<vector<int>> cond1;
	vector<vector<int>> cond2;

	vector<double> local_cond;
	vector<double> local_cond_vl;

	//матрицы и вектора
	vector<vector<double>> G;//локальные матрицы в плотном формате
	vector<vector<double>> M;
	vector<vector<double>> C;
	vector<vector<double>> p_loc, c_loc;//для блочного
	vector<vector<double>> local;
	vector<double> b_loc;
	vector<double> fc;//локальный вектор
	vector<double> fs;//локальный вектор
	vector<double> f;//вспомогательный, значения функции в координатах

	vector<vector<double>> G_small, M_small;

	vector<double> x_loc, y_loc, z_loc;//координаты локальных узлов каждого элемента

	//глобальные матрицы в разреженном формате
	vector<double> ggu, ggl;//элементы
	vector<double>ggu_new, ggl_new;
	vector<double> di;//диагональные
	vector<double> di_new;
	vector<int> ig;//портрет
	vector<int> jg;//портрет

	//матрица в профильном
	vector<int> ig_pr;
	vector<double> di_pr;
	vector<double> ggu_pr, ggl_pr;

	//глобальные вектора
	vector<double> b; //глобальная пр часть
	vector<double> q; //истинное решение

	//slae
	vector<double> r, p, s, temp1, temp2, r0, t, y_temp, temp3;

	double BSG_alfa, BSG_beta, BSG_gamma;
	double st_time, end_time, prog_time;

	int flag;
	//vector<int> global_num;


	SLAU(string net_file, string cond_file_1, int _flag) {
		flag = _flag;
		read_net(net_file);
		init_net();
		read_cond(cond_file_1);
		resize_func();
	}

	~SLAU() {

		for (int i = 0; i < 8; i++) {
			G[i].~vector();
			M[i].~vector();
			C[i].~vector();
			local[i].~vector();
			p_loc[i].~vector();
			c_loc[i].~vector();
		}

		G.~vector();
		M.~vector();
		C.~vector();
		local.~vector();
		p_loc.~vector();
		c_loc.~vector();

		b_loc.~vector();
		fs.~vector();
		fc.~vector();
		f.~vector();
		//global_num.resize(2);

		ig.~vector();

		b.~vector();
		q.~vector();

		x_loc.~vector();
		y_loc.~vector();
		z_loc.~vector();

		for (int i = 0; i < 2; i++) {
			G_small[i].~vector();
			M_small[i].~vector();
		}
		G_small.~vector();
		M_small.~vector();

		//для Слау
		r.~vector();
		z.~vector();
		p.~vector();
		temp1.~vector();
		temp2.~vector();
		if (flag == 2) {
			r0.~vector();
			temp3.~vector();
			t.~vector();
			y_temp.~vector();
		}
	}

	void read_net(string net_file);
	void init_net();
	void resize_func();
	void read_cond(string cond_file_1);
	void write_small_matrix();
	int get_global(int r, int s, int p, int j);
	double get_Ugs(double x, double y, double z);
	double get_Ugc(double x, double y, double z);
	double func_s(int num_area, double x, double y, double z);
	double func_c(int num_area, double x, double y, double z);
	void make_portrait();
	void make_local_coor(int r, int s, int p);
	void make_local_matrix(double hx, double hy, double hz);
	void make_local_vec(int num_area);
	void clean_matrix();
	void make_global();
	void add_first();
	void add_local(int r, int s, int p);

	//SLAE 
	void LU_sq();
	void mult(vector<double>& MV, vector<double> vec);
	double skal_mult(vector<double> vec1, vector<double> vec2);
	long double norm(vector<double> vec);
	void mult_pr(vector<double> aa, vector<double>di, vector<double>& y, vector<double> bb);
	void mult_obr(vector<double> aa, vector<double>di, vector<double>& y, vector<double> bb);
	void multU(vector<double>& MV, vector<double>vec);
	void mult_tr(vector<double>& MV, vector<double> vec);
	void LOS_sq();
	void generate_profil();
	void LU_mult_pr();
	void LU_mult_obr();
	void LU();
	void LU_sol();
	void BSG_Stab_LU_sq();
	void do_smth();
};