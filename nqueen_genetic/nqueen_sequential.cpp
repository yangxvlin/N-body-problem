/*
Alex Marroquin
CSCI 4350 Artificial Intelligence
Homework 3: N-Queen Solution with Genetic Algorithm
main.cpp
*/

#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<time.h>
#include<vector>

using namespace std;

const int POP_SIZE = 100;
int fit_p1, fit_p2;
ofstream out("output.txt");

//Checks if number already exists in vector
bool checkExists(vector<int> a, int spot, int val)
{
	for(int i = 0; i < spot; i++)
	{
		if(a[i] == val)
			return true;
	}
	return false;
}
//Calculates the number of collisions
int numCollision(vector<int> a, int num)
{
	int collisionCount = 0;
	bool upCollision = false;
	bool downCollision = false;
	for(int i = 0; i < num; i++)
	{
		for(int j = i+1; j < num; j++)
		{
			if(j < num)
			{
				int dx = i - j;
				int dy = a[i] - a[j];
				if(abs(dx) == abs(dy) && (!upCollision || !downCollision))
				{
					if(dy > 0)
						upCollision = true;
					else
						downCollision = true;
					collisionCount++;
				}
			}
		}
		upCollision = false;
		downCollision = false;
	}
	return collisionCount;
}

//Used for storing the piece arrangments and fitness of
//that particular arrangement
class Queen
{
public:
	double fitness;
	vector<int> places;
	Queen();
	Queen(int);
	void print();
	void outPrint();
};
//Outputs queen arrangement data
void Queen::print()
{
	cout<<endl;
	for(int i = 0; i < places.size(); i++)
	{
		cout<<places[i]<<" ";
	}
	cout<<endl;
}
//Outputs to file queen arrangement data
void Queen::outPrint()
{
	for(int i = 0; i < places.size(); i++)
	{
		out<<places[i]<<" ";
	}
}
Queen::Queen()
{
	/*for(int i = 0; i < 10; i++)
	{
		places.push_back(i);
	}*/
}
Queen::Queen(int num)
{
	//Initiates positions
	int r; 
	//int collisionCount = 0;
	for(int i = 0; i < num; i++)
	{
		do{
			r = rand()%num + 1;
		}while(checkExists(places,i,r));

		places.push_back(r);
	}
	//print();

	//Calculates the fitness of the arrangement
	int collisionCount = numCollision(places,num);
	fitness = 1/(collisionCount + 0.01);
}

//Stores an entire populatin of arrangements
class Population
{
public:
	Queen *population[POP_SIZE];
	Population();
	Population(int);
};
Population::Population()
{
	for(int i = 0; i < 10; i++)
	{
		population[i] = new Queen(8);
		population[i]->print();
	}
}
Population::Population(int size)
{
	for(int i = 0; i < POP_SIZE; i++)
	{
		population[i] = new Queen(size);
	}
}

//Find most fit in population
Queen* findMostFit(Population *pop, int size)
{
	Queen *a = pop->population[0];
	for(int i = 0; i < size; i++)
	{
		if(a != NULL && pop->population[i] != NULL && a->fitness < pop->population[i]->fitness)
		{
			a = pop->population[i];
			fit_p1 = i;
		}
	}

	return a;
}
Queen* findMostFit(vector<Queen*> pop, int size)
{
	Queen *a = pop[0];
	for(int i = 0; i < size; i++)
	{
		if(a != NULL && pop[i] != NULL && a->fitness < pop[i]->fitness)
		{
			a = pop[i];
			fit_p2 = i;
		}
	}

	return a;
}

//Get position of most fit
int mostFitPos(Population *pop, int size)
{
	int pos = 0;
	for(int i = 0; i < size; i ++)
	{
		if(pop->population[pos] < pop->population[i])
			pos = i;
	}
	return pos;
}
int mostFitPos(vector<Queen*> pop, int size)
{
	int pos = 0;
	for(int i = 0; i < size; i ++)
	{
		if(pop[pos] < pop[i])
			pos = i;
	}
	return pos;
}

//Method for selecting parents, with substitution
bool getMate(vector<Queen*> *a, Population h)
{
	int r1, r2, r3;
	int current_member = 1;
	Queen *b;
	while(current_member <= 40)
	{
		//Select random members
		r1 = rand()%POP_SIZE;
		do{
			r2 = rand()%POP_SIZE;
		}while(r2 == r1);
		do{
			r3 = rand()%POP_SIZE;
		}while(r3 == r1 || r3 == r2);

		//Find most fit among population
		if(h.population[r1]->fitness > h.population[r2]->fitness && h.population[r1]->fitness > h.population[r3]->fitness)
			b = h.population[r1];
		else if(h.population[r2]->fitness > h.population[r1]->fitness && h.population[r2]->fitness > h.population[r3]->fitness)
			b = h.population[r2];
		else
			b = h.population[r3];
		
		//Add best to mating pool and increase counter
		a->push_back(b);
		current_member++;
	}
	return true;

}

int main()
{
	int quant; //Desired N value
	cout<<"Welcome to the N-Queen solver. How much do you want N to be?  ";
	cin>>quant;
	bool found = false; //If solution is found
	int generations = 0; //How many generations til solution

	out<<"Runs"<<"		"<<"Generations"<<"	"<<"Optimal Fitness"<<"		"<<"Solution Vector"<<endl;

	//Initialize population
	srand (time(NULL));
	system("CLS");

	for(int times = 0; times < 50; times++)
	{
		Population herd(quant); //Population which algorithm is applied to
		while(!found)
		{
			//Parent Selection
			vector<Queen*> mategroup; //Mating group
			getMate(&mategroup,herd);

			//Reproduction
			/***************************
			Copies half of one parent to one child, then copies rest of parent
			***************************/
			vector<Queen*> children;
			for(int i = 0; i < 20; i++)
			{
				//Select parents
				Queen *p1 = mategroup[i];
				Queen *p2 = mategroup[i+1];

				//Initialize Children
				Queen *c1 = new Queen;
				Queen *c2 = new Queen;

				//Add half of one parent to one child
				for(int j = 0; j < quant/2; j++)
				{
					c1->places.push_back(p1->places[j]);
					c2->places.push_back(p2->places[j]);
				}

				//Add rest of positions
				int cnt1 = quant/2;
				int cnt2 = cnt1;
				for(int k = 0; k < quant; k++)
				{
					if(!checkExists(c1->places,cnt1,p2->places[k]))
					{
						c1->places.push_back(p2->places[k]);
						cnt1++;
					}
					if(!checkExists(c2->places,cnt2,p1->places[k]))
					{
						c2->places.push_back(p1->places[k]);
						cnt2++;
					}
				}
		
				//Check children values
				/*for(int k = 0; k < quant; k++)
				{
					cout<<c1->places[k]<<" ";
				}
				cout<<"		";
				for(int k = 0; k < quant; k++)
				{
					cout<<c2->places[k]<<" ";
				}
				cout<<endl<<endl;*/

				//Mutation
				/**********************************************************
				Select a random number between 0 and 3.  If the number is 1
				then no mutation occurs, else two spots at random are selected
				and switched
				**********************************************************/
				int mut = rand()%4;
				if(mut != 1)
				{
					int pos1 = rand()%quant;
					int pos2;
					do{
						pos2 = rand()%quant;
					}while(pos1 == pos2);

					int val1 = c1->places[pos1];
					int val2 = c1->places[pos2];

					c1->places[pos1] = val2;
					c1->places[pos2] = val1;
				}
				mut = rand()%4;
				if(mut != 1)
				{
					int pos1 = rand()%quant;
					int pos2;
					do{
						pos2 = rand()%quant;
					}while(pos1 == pos2);

					int val1 = c2->places[pos1];
					int val2 = c2->places[pos2];

					c2->places[pos1] = val2;
					c2->places[pos2] = val1;
				}

				//Check children values after mutation
				/*for(int k = 0; k < quant; k++)
				{
					cout<<c1->places[k]<<" ";
				}
				cout<<"		";
				for(int k = 0; k < quant; k++)
				{
					cout<<c2->places[k]<<" ";
				}*/

				//Calculate fitness of children
				int col1 = numCollision(c1->places,quant);
				int col2 = numCollision(c2->places,quant);
				c1->fitness = 1/(col1 + 0.01);
				c2->fitness = 1/(col2 + 0.01);
				
				if(c1->fitness == 100 || c2->fitness == 100)
					found = true;

				children.push_back(c1);
				children.push_back(c2);
			}
			generations++;

			system("CLS");

			//Survivor Selection
			/**********************************************************
			The 200 best fit individuals are found and stored in the dummy variable
			best.  Afterwards they are copied to herd.
			***********************************************************/
			Population best(quant); //Temporarily stores best individuals
			Queen *fit1, *fit2; //To store best in population and children, respectively
			//int fit_p1, fit_p2;
			int v_cnt = 0;
			int h_cnt = 0;
			for(int i = 0; i < POP_SIZE; i++)
			{
				fit1 = findMostFit(&herd,POP_SIZE - i);
				fit2 = findMostFit(children,40);

				//Find which is most fit
				if(fit2 == NULL || fit1->fitness > fit2->fitness)
				{
					best.population[i] = fit1;
					cout<<"Selected fit1\n";

					herd.population[fit_p1] = herd.population[POP_SIZE - 1 - h_cnt];
					herd.population[POP_SIZE - 1 - h_cnt] = NULL;
					h_cnt++;
				}
				else
				{
					best.population[i] = fit2;
					cout<<"Selected fit2\n";
					children[fit_p2] = children[39 - v_cnt];
					children[39 - v_cnt] = NULL;
					v_cnt++;
				}
				system("CLS");
			}
			herd = best;
			cout<<"Finished cycle\n";
			system("CLS");
		}//end while loop
		out<<times+1<<"		"<<generations<<"		"<<herd.population[0]->fitness<<"			";
		herd.population[0]->outPrint();
		out<<endl;

		system("CLS");
		found = false;
		generations = 0;

	}//end for loop

	out.close();
	system("pause");
	return 0;
}

// g++ -O3 nqueen_sequential.cpp -o nqueen
// ./nqueen