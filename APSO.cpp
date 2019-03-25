#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>
#define pi acos(-1)
#define E  2.71828182845904523536

using namespace std;

const int Exploration = 0;
const int Exploitation = 1;
const int Convergence = 2;
const int Jumping_out = 3;

double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}
class Particle
{
public:
	vector<double> position;
	vector<double>velocity;
	vector<double>pBest;

	Particle() {}

	Particle(vector<double> position, vector<double>velocity, vector<double>best_position, double best_fitness)
	{
		this->position = position;
		this->velocity = velocity;
		this->pBest = best_position;
	}
};

bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}

class PSO
{
public:

	PSO(int dim, int m, int Tmax, double max, double min, double c1, double c2, double wmax, double wmin, double dt, double percent, 
		double Varmax, double Varmin)
	{
		this->dim = dim;
		this->m = m;
		this->Tmax = Tmax;
		this->max = max;
		this->min = min;
		this->c1 = c1;
		this->c2 = c2;
		this->wmax = wmax;
		this->wmin = wmin;
		this->dt = dt;
		this->percent = percent;
		this->Varmax = Varmax;
		this->Varmin = Varmin;
		number_increase_c1 = 0;
		number_increase_c2 = 0;
		number_s1 = 0;
		number_s2 = 0;
		number_s3 = 0;
		number_s4 = 0;
		w = 0.9;
		state = Exploration;
		particles.resize(m);
		gBestIndex = 0;
		d.resize(m);
	}

	void setGbest()
	{
		for (int i = 0; i < m; i++)
		{
			if (fitnessFunction(particles[i].pBest) < fitnessFunction(particles[gBestIndex].pBest))
			{
				gBestIndex = i;
			}
		}

	}

	double fitnessFunction(vector<double>pos)
	{
		double result = 0.0;

		for (int i = 0; i < dim; i++)
		{
			result += pow(pos[i], 2) + 10 - 10*cos(2*pi*pos[i]);
		}
		return result;
	}

	void initialParticles(int i)
	{
		particles[i].position.resize(dim);
		particles[i].velocity.resize(dim);
		particles[i].pBest.resize(dim);
		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);
			particles[i].position[j] = randDouble(min, max);
			particles[i].velocity[j] = randDouble(-range, range);
			particles[i].pBest[j] = particles[i].position[j];
		}
	}

	void initialAllParticles()
	{
		for (int i = 0; i < m; i++)
		{
			initialParticles(i);
		}
		setGbest();
	}

	void caculateDistance(int i)
	{
		double dis = 0.0;
		double temp;
		for (int j = 0; j < m; j++)
		{
			temp = 0.0;
			for (int k = 0; k < dim; k++)
			{
				temp += pow((particles[i].position[k] - particles[j].position[k]), 2);
			}
			dis += sqrt(temp);
		}
		d[i] = dis / (m - 1);
		
	}

	void caculateAllDistance()
	{
		for (int i = 0; i < m; i++)
		{
			caculateDistance(i);
		}

		int max_index = 0;
		int min_index = 0;
		for (int i = 0; i < m; i++)
		{
			if (d[i] < d[min_index])
			{
				min_index = i;
			}
			else if (d[i] > d[max_index])
			{
				max_index = i;
			}
		}
		dmax = d[max_index];
		dmin = d[min_index];

		dg = d[gBestIndex];
	}

	void evolutionaryFactor()
	{
		f = (dg - dmin) / (dmax - dmin);
		//cout << f << endl;
		//result.push_back(f);
	}

	void inertiaWeight()
	{
		w = 1 / (1 + 1.5*exp(-2.6*f));
		result1.push_back(w);
		result2.push_back(f);
	}

	double f1()
	{
		if (f >= 0 && f <= 0.4)
			return 0;
		if (f > 0.4 && f <= 0.6)
			return (5*f - 2);
		if (f > 0.6 && f <= 0.7)
		{
			state = Exploration;
			return 1;
		}
		if (f > 0.7 && f <= 0.8)
			return (-10*f + 8);
		if (f > 0.8 && f <= 1)
			return 0;
	}

	double f2()
	{
		if (f >= 0 && f <= 0.2)
			return 0;
		if (f > 0.2 && f <= 0.3)
			return (10 * f - 2);
		if (f > 0.3 && f <= 0.4)
		{
			state = Exploitation;
			return 1;
		}
		if (f > 0.4 && f <= 0.6)
			return (-5 * f + 3);
		if (f > 0.6 && f <= 1)
			return 0;
	}

	double f3()
	{
		if (f >= 0 && f <= 0.1)
		{
			state = Convergence;
			return 1;
		}
		if (f > 0.1 && f <= 0.3)
			return (5 * f - 1.5);
		
		if (f > 0.3 && f <= 1)
			return 0;
	}

	double f4()
	{
		if (f >= 0 && f <= 0.7)
			return 0;
		if (f > 0.7 && f <= 0.9)
			return (5 * f - 3.5);
		if (f > 0.9 && f <= 1)
		{
			state = Jumping_out;
			return 1;
		}
	}

	void identifyState()
	{
		vector<double>st;
		bool end = false;
		st.resize(4);
		st[0] = f1();
		st[1] = f2();
		st[2] = f3();
		st[3] = f4();
		/*
		int ma = 0;
		for (int i = 0; i < 4; i++)
		{
			if (st[i] > st[ma])
				ma = i;
		}
		state = ma;
		//cout << "state " << state + 1 << endl;
		*/
		int number_of_zero = 0;
		for (int i = 0; i < st.size(); i++)
		{
			if (st[i] == 1)
			{
				cout << f << endl;
				cout << state << endl;
				return;
			}
			if (st[i] == 0)
				number_of_zero++;
		}

		int s;
		int pre_state = state;

		if (number_of_zero == 3)
		{
			for (int i = 0; i < st.size(); i++)
			{
				if (st[i] != 0)
					s = i;
			}
			state = s;
		}
		else if(number_of_zero == 2)
		{
			int s1 = -1, s2 = -1;
			for (int i = 0; i < st.size(); i++)
			{
				if (st[i] != 0 && s1 == -1)
				{
					s1 = i;
				}
				else if (st[i] != 0 && s1 != -1)
				{
					s2 = i;
				}
			}

			
			if (s1 == Exploration && s2 == Exploitation)
			{
				if (pre_state == Exploration)
				{
					state = Exploration;
				}
				else if (pre_state == Jumping_out)
				{
					state = Exploration;
				}
				else
				{
					state = Exploitation;
				}
			}

			if (s1 == Exploration && s2 == Jumping_out)
			{
				if (pre_state == Exploration)
				{
					state = Exploration;
				}
				else if (pre_state == Exploitation)
				{
					state = Exploration;
				}
				else
				{
					state = Jumping_out;
				}
			}

			else if (s1 == Exploitation && s2 == Convergence)
			{
				if (pre_state == Exploitation)
				{
					state = Exploitation;
				}
				else if (pre_state == Exploration)
				{
					state = Exploitation;
				}
				else
				{
					state = Convergence;
				}
			}
		}
		cout << f << endl;
		cout << state << endl;
		
	}

	void Elite_Learning()
	{
		static default_random_engine engine(time(nullptr));
		uniform_int_distribution<int> dis(0, dim-1);
		int d = dis(engine);
		Var = Varmax - (Varmax - Varmin)*T / Tmax;

		normal_distribution<double>gaussin(0, Var);
		double n = gaussin(engine);

		Particle p = particles[gBestIndex];
		p.pBest[d] = p.pBest[d] + (max - min)*n;
		if (p.pBest[d] > max)
		{
			p.pBest[d] = max;
		}
		else if (p.pBest[d] < min)
		{
			p.pBest[d] = min;
		}

		if (fitnessFunction(p.pBest) <= fitnessFunction(particles[gBestIndex].pBest))
		{
			particles[gBestIndex] = p;
		}
		else
		{
			int index = 0;
			for (int i = 1; i < particles.size(); i++)
			{
				if (fitnessFunction(particles[index].pBest) < fitnessFunction(particles[i].pBest))
					index = i;
			}
			//particles[index].pBest = p.pBest;

			if (fitnessFunction(p.pBest) < fitnessFunction(particles[index].pBest))
				particles[index] = p;
			
		}
		

	}


	void updateC()
	{
		double pre_c1 = c1;
		double pre_c2 = c2;
		double n1 = randDouble(0.05, 0.1);
		double n2 = randDouble(0.05, 0.1);
		double d1 = randDouble(0, n1);
		double d2 = randDouble(0, n2);
		
		if (state == Exploration)
		{
			c1 = c1 + d1;
			c2 = c2 - d2;
			number_s1++;
		}

		else if (state == Exploitation)
		{
			c1 = c1 + 0.5*d1;
			c2 = c2 - 0.5*d2;
		}

		else if (state == Convergence)
		{
			c1 = c1 + 0.5*d1;
			c2 = c2 + 0.5*d2;
			Elite_Learning();
		}

		else if (state == Jumping_out)
		{
			c1 = c1 - d1;
			c2 = c2 + d2;
		}

		if (c1 > 2.5)
			c1 = 2.5;
		else if (c1 < 1.5)
			c1 = 1.5;
		if (c2 > 2.5)
			c2 = 2.5;
		else if (c2 < 1.5)
			c2 = 1.5;

		double temp1 = c1;
		double temp2 = c2;
		if (c1 + c2 > 4.0)
		{
			
			c1 = temp1 / (temp1 + temp2)*4.0;
			c2 = temp2 / (temp1 + temp2)*4.0;
		}

		result3.push_back(c1);
		result4.push_back(c2);
	}

	void caculatePSD()
	{
		vector<double>mean_position;
		mean_position.resize(dim);
		
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < m; j++)
			{
				mean_position[i] += particles[j].position[i];
			}
			mean_position[i] = mean_position[i] / m;
			//cout << mean_position[i] << endl;
		}

		double temp = 0.0;
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				temp += pow((particles[i].position[j] - mean_position[j]), 2);
			}
			
		}
		temp = sqrt(temp / (m - 1));

		result5.push_back(temp);

		result6.push_back(fitnessFunction(particles[gBestIndex].pBest));
	}

	void updateParticle(int i)
	{
		for (int j = 0; j < dim; j++)
		{
			double last_position = particles[i].position[j];
			double range = percent * (max - min);

			

			particles[i].velocity[j] = w * particles[i].velocity[j] +
				c1 * randDouble(0, 1) * (particles[i].pBest[j] - particles[i].position[j])
				+ c2 * randDouble(0, 1) * (particles[gBestIndex].pBest[j] - particles[i].position[j]);

			if (particles[i].velocity[j] > range)
				particles[i].velocity[j] = range;

			else if (particles[i].velocity[j] < -range)
				particles[i].velocity[j] = -range;


			particles[i].position[j] += dt * particles[i].velocity[j];

			if (particles[i].position[j] > max)
			{
				particles[i].position[j] = max;
			}
			else if (particles[i].position[j] < min)
			{
				particles[i].position[j] = min;
			}
			
		}

		if (fitnessFunction(particles[i].position) < fitnessFunction(particles[i].pBest))
		{
			particles[i].pBest = particles[i].position;
		}

	}


	void updateAllParticles()
	{
		for (int i = 0; i < m; i++)
		{
			updateParticle(i);
		}
		setGbest();
		T++;
	}

	double getFitness(int i)
	{
		return fitnessFunction(particles[i].pBest);
	}
	double getFitness()
	{
		return fitnessFunction(particles[gBestIndex].pBest);
	}

	int number_increase_c1;
	int number_increase_c2;
	int number_s1;
	int number_s2;
	int number_s3;
	int number_s4;

	vector<double>result1;
	vector<double>result2;
	vector<double>result3;
	vector<double>result4;
	vector<double>result5;
	vector<double>result6;

private:
	int dim;
	int m;//number of instances

	int T;
	int Tmax;

	double w;
	double max;
	double min;
	double c1;
	double c2;
	double wmax;
	double wmin;
	double Varmax;
	double Varmin;
	double Var;

	double dt;//时间步长
	double percent;

	int gBestIndex;

	vector<Particle> particles;

	double f;
	vector<double>d;
	double dmax;
	double dmin;
	double dg;

	int state;

	
	
};

void run(vector<double>& result1, vector<double>& result2, vector<double>& result3,
	vector<double>& result4, vector<double>& result5, vector<double>& result6)
{
	int dim = 30;
	int m = 20;
	int Tmax = 2000;
	double max = 5.12;
	double min = -5.12;
	double c1 = 2;
	double c2 = 2;
	double wmax = 0.9;
	double wmin = 0.4;
	double dt = 1.0;
	double percent = 0.2;
	double Varmax = 1.0;
	double Varmin = 0.1;

	PSO pso = PSO(dim, m, Tmax, max, min, c1, c2, wmax, wmin, dt, percent, Varmax, Varmin);
	pso.initialAllParticles();


	vector<double>fitness;
	fitness.push_back(pso.getFitness());

	for (int i = 0; i < Tmax; i++)
	{
		pso.caculateAllDistance();
		pso.evolutionaryFactor();
		pso.identifyState();
		pso.updateC();

		pso.inertiaWeight();
		pso.updateAllParticles();
		pso.caculatePSD();
		//fitness.push_back(pso.getFitness());
		if ((1 + i) % 100 == 0) {
			fitness.push_back(pso.getFitness());

		}
		cout << "第" << i << "次迭代结果：";
		cout << ", fitness = " << pso.getFitness() << endl;



	}
	result1 = pso.result1;
	result2 = pso.result2;
	result3 = pso.result3;
	result4 = pso.result4;
	result5 = pso.result5;
	result6 = pso.result6;
}

int main()
{
	vector<double> result1;
	vector<double> result2;
	vector<double> result3;
	vector<double> result4;
	vector<double> result5;
	vector<double> result6;

	
	vector<double> result1_temp;
	vector<double> result2_temp;
	vector<double> result3_temp;
	vector<double> result4_temp;
	vector<double> result5_temp;
	vector<double> result6_temp;

	int test_time = 5;
	int gape = 5;

	run(result1, result2, result3, result4, result5, result6);
	for (int i = 1; i < test_time; i++)
	{
		run(result1_temp, result2_temp, result3_temp, result4_temp, result5_temp, result6_temp);
		for (int j = 0; j < result1.size(); j++)
		{
			result1[j] += result1_temp[j];
			result2[j] += result2_temp[j];
			result3[j] += result3_temp[j];
			result4[j] += result4_temp[j];
			result5[j] += result5_temp[j];
			result6[j] += result6_temp[j];
		}
	}

	for (int j = 0; j < result1.size(); j++)
	{
		result1[j] = result1[j] / test_time;
		result2[j] = result2[j] / test_time;
		result3[j] = result3[j] / test_time;
		result4[j] = result4[j] / test_time;
		result5[j] = result5[j] / test_time;
		result6[j] = result6[j] / test_time;
	}
	cout << endl;
	for (int j = 0; j < result1.size(); j++)
	{
		if(j%gape ==0)
			cout << result1[j] << " ";
	}
	cout << endl;
	for (int j = 0; j < result1.size(); j++)
	{
		if (j % gape == 0)
		cout << result2[j] << " ";
	}
	cout << endl;
	for (int j = 0; j < result1.size(); j++)
	{
		if (j % gape == 0)
		cout << result3[j] << " ";
	}
	cout << endl;
	for (int j = 0; j < result1.size(); j++)
	{
		if (j % gape == 0)
		cout << result4[j] << " ";
	}
	cout << endl;
	for (int j = 0; j < result1.size(); j++)
	{
		if (j % gape == 0)
		cout << result5[j] << " ";
	}
	cout << endl;
	for (int j = 0; j < result1.size(); j++)
	{
		if (j % gape == 0)
		cout << result6[j] << " ";
	}
	cout << endl;
	
	system("pause");
}