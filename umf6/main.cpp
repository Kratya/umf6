#include "slau.h"

int main() {


	ofstream res("LU10-4.txt");
	SLAU slv = SLAU("net", "cond1", 2);


	cout << slv.ggl.max_size();
	cout << slv.ggu_new.max_size();


	slv.omega = 1e-4;

	res << "===========lyam=1e+2============" << endl;
	cout << "===========lyam=1e+2============" << endl;
	slv.lyam = 100;
	{
		cout << "__________sigma=0_____________" << endl;
		res << "__________sigma=0_____________" << endl;
		slv.sigma = 0;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl;
	}
	{
		cout << "__________sigma=1e4_____________" << endl;
		res << "__________sigma=1e4_____________" << endl;
		slv.sigma = 1e+4;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl; }
	{
		cout << "__________sigma=1e8_____________" << endl;
		res << "__________sigma=1e8_____________" << endl;
		slv.sigma = 1e+8;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl; }

	cout << "===========lyam=1e+4=========" << endl;
	res << "===========lyam=1e+4=========" << endl;
	slv.lyam = 1e+4;
	{
		cout << "__________sigma=0_____________" << endl;
		res << "__________sigma=0_____________" << endl;
		slv.sigma = 0;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl;
	}
	{
		cout << "__________sigma=1e4_____________" << endl;
		res << "__________sigma=1e4_____________" << endl;
		slv.sigma = 1e+4;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl; }
	{
		cout << "__________sigma=1e8_____________" << endl;
		res << "__________sigma=1e8_____________" << endl;
		slv.sigma = 1e+8;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl; }

	cout << "============lyam=8e+5==============" << endl;
	res << "============lyam=8e+5==============" << endl;
	slv.lyam = 8e+5;
	{
		cout << "__________sigma=0_____________" << endl;
		res << "__________sigma=0_____________" << endl;
		slv.sigma = 0;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl;
	}
	{
		cout << "__________sigma=1e4_____________" << endl;
		res << "__________sigma=1e4_____________" << endl;
		slv.sigma = 1e+4;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl; }
	{
		cout << "__________sigma=1e8_____________" << endl;
		res << "__________sigma=1e8_____________" << endl;
		slv.sigma = 1e+8;
		res << endl << "1e-11" << endl;
		slv.hi = 1e-11;
		slv.do_smth();

		res << "time=" << slv.prog_time / 1000 << endl << endl;
		res << "nev=" << slv.nev << endl;
		res << endl << "iter=" << slv.iter << endl; }

	res.close();

	
	ofstream res2("BCG10+5.txt");
	SLAU slv2 = SLAU("net", "cond1", 2);

	slv2.omega = 1e+5;

	res2 << "===========lyam=1e+2============" << endl;
	cout << "===========lyam=1e+2============" << endl;
	slv2.lyam = 100;
	{
		cout << "__________sigma=0_____________" << endl;
		res2 << "__________sigma=0_____________" << endl;
		slv2.sigma = 0;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl;
	}
	{
		cout << "__________sigma=1e4_____________" << endl;
		res2 << "__________sigma=1e4_____________" << endl;
		slv2.sigma = 1e+4;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl; }
	{
		cout << "__________sigma=1e8_____________" << endl;
		res2 << "__________sigma=1e8_____________" << endl;
		slv2.sigma = 1e+8;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl; }

	cout << "===========lyam=1e+4=========" << endl;
	res2 << "===========lyam=1e+4=========" << endl;
	slv2.lyam = 1e+4;
	{
		cout << "__________sigma=0_____________" << endl;
		res2 << "__________sigma=0_____________" << endl;
		slv2.sigma = 0;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl;
	}
	{
		cout << "__________sigma=1e4_____________" << endl;
		res2 << "__________sigma=1e4_____________" << endl;
		slv2.sigma = 1e+4;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl; }
	{
		cout << "__________sigma=1e8_____________" << endl;
		res2 << "__________sigma=1e8_____________" << endl;
		slv2.sigma = 1e+8;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl; }

	cout << "============lyam=8e+5==============" << endl;
	res2 << "============lyam=8e+5==============" << endl;
	slv2.lyam = 8e+5;
	{
		cout << "__________sigma=0_____________" << endl;
		res2 << "__________sigma=0_____________" << endl;
		slv2.sigma = 0;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl;
	}
	{
		cout << "__________sigma=1e4_____________" << endl;
		res2 << "__________sigma=1e4_____________" << endl;
		slv2.sigma = 1e+4;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl; }
	{
		cout << "__________sigma=1e8_____________" << endl;
		res2 << "__________sigma=1e8_____________" << endl;
		slv2.sigma = 1e+8;
		res2 << endl << "1e-11" << endl;
		slv2.hi = 1e-11;
		slv2.do_smth();

		res2 << "time=" << slv2.prog_time / 1000 << endl << endl;
		res2 << "nev=" << slv2.nev << endl;
		res2 << endl << "iter=" << slv2.iter << endl; }

	res2.close();



	ofstream fout("out.txt");
for (int i = 0; i < slv.m * 2; i++)
	fout << slv.q[i] << endl;
	res << "time=" << slv.end_time/1000 ;

	system("pause");
	return 0;

}