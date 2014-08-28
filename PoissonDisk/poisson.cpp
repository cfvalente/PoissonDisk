#include <iostream>
#include <vector>
#include <list>
#include <random>
#include <time.h>
#include <omp.h>

#define PI 3.14159265359

using namespace std;

struct point
{
	double x;
	double y;
	double z;
};


default_random_engine generator;


double cell_size;
int ***grid;

vector<int> rlist;
vector<struct point> points;

int rand_int(int vi, int vf)
{
	uniform_int_distribution<int> distribution(vi,vf);
	return distribution(generator);
}

double rand_double(int vi, int vf)
{

	uniform_real_distribution<double> distribution(vi,vf);
	return distribution(generator);
}

void new_grid(int width, int height, int length)
{
	int grid_w,grid_h,grid_l;
	grid_w = ceil(width/cell_size);
	grid_h = ceil(height/cell_size);
	grid_l = ceil(length/cell_size);
	grid = new int** [grid_w];
#pragma omp parallel for schedule(static)
	for(int i=0; i<=grid_w; i++)
	{
		grid[i] = new int* [grid_h];
		for(int j=0; j<=grid_h; j++)
		{
			grid[i][j] = new int [grid_l];
			for(int k=0; k<=grid_l; k++)
			{
				grid[i][j][k] = -1;
			}
		}
	}
}

void random_list_insert(int x)
{
	rlist.push_back(x);
}

int random_list_size()
{
	return rlist.size();
}

int random_list_get()
{
	if(rlist.size() > 0)
	{
		int rand = rand_int(0,rlist.size()-1);
		int element = rlist[rand];
		rlist.erase(rlist.begin()+rand);
		return element;
	}
	return -1;
}

bool random_list_empty()
{
	if(rlist.size() > 0) return 0;
	return 1;
}

void insert_in_grid(struct point p, int ind)
{
	int x,y,z;
	x = floor(p.x/cell_size);
	y = floor(p.y/cell_size);
	z = floor(p.z/cell_size);
	grid[x][y][z] = ind;
}

struct point create_random_point(struct point p, double min_dist)
{
	struct point generated_point;
	double r1 = rand_double(0,1); 
	double r2 = rand_double(0,1); 
	double r3 = rand_double(0,1); 
	double radius = min_dist * (r1 + 1);
	double theta = 2 * PI * r2;
	double phi = 2 * PI * r3;

	generated_point.x = p.x+radius*cos(theta)*sin(phi);
	generated_point.y = p.y+radius*sin(theta)*sin(phi);
	generated_point.z = p.z+radius*cos(phi);

	return generated_point;
}

bool valid_point(struct point p, int width, int heigth, int length)
{
	if(p.x >= 0 && p.x <= width)
		if(p.y >= 0 && p.y <= heigth)
			if(p.z >= 0 && p.z <= length)
				return 1;
	return 0;
}

double point_distance(struct point p1, struct point p2)
{
	return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
}

bool no_close_neighbours(struct point new_point, int width, int heigth, int length, double min_dist)
{
	int npx_grid,npy_grid,npz_grid;
	bool close_neighbour = 0;
	npx_grid = floor(new_point.x/cell_size);
	npy_grid = floor(new_point.y/cell_size);
	npz_grid = floor(new_point.z/cell_size);
#pragma omp parallel for schedule(static, 1)
	for(int x = npx_grid-2; x <= npx_grid+2; x++)
	{
		if(x >=0 && x <= ceil(width/cell_size) && !close_neighbour)
		{
			for(int y = npy_grid-2; y <= npy_grid+2; y++)
			{
				if(y >=0 && y <= ceil(heigth/cell_size) && !close_neighbour)
				{
					for(int z = npz_grid-2; z <= npz_grid+2; z++)
					{
						if(z >=0 && z <= ceil(length/cell_size) && !close_neighbour)
						{
							if(grid[x][y][z]!=-1)
							{
								double dist;
								dist = point_distance(new_point, points[grid[x][y][z]]);
								if(dist < min_dist)
									close_neighbour = 1;
							}
						}
					}
				}
			}
		}
	}
	if(!close_neighbour)
		return 1;
	return 0;
}

void init(int width, int heigth, int length, double min_dist)
{
	struct point p;
	int ind;
	cell_size = min_dist/sqrt(2.0);
	new_grid(width, heigth,length);
	generator.seed(time(NULL));

	p.x = rand_double(0, width);
	p.y = rand_double(0, heigth);
	p.z = rand_double(0, length);
	points.push_back(p);
	ind = points.size()-1;
	random_list_insert(ind);
	insert_in_grid(p,ind);
}

void process_list(int width, int heigth, int length, double min_dist, int point_count)
{
	struct point p,new_point;
	while(!random_list_empty())
	{
		p = points[random_list_get()];
		for(int i = 0; i < point_count; i++)
		{
			new_point = create_random_point(p, min_dist);
			if(valid_point(new_point, width, heigth, length))
			{
				if(no_close_neighbours(new_point, width, heigth, length, min_dist))
				{
					int ind;
					points.push_back(new_point);
					ind = points.size()-1;
					random_list_insert(ind);
					insert_in_grid(new_point,ind);
				}
			}
		}

	}
}

void create_spheres_pbrt()
{
	for(int  i = 0; i<points.size(); i++)
	{
		cout << "P[" << i << "] -- X: " << points[i].x << "  - Y: " << points[i].y << "  - Z: " << points[i].z << endl;
	}
}

int main(int argc, char *argv[])
{
	int width, heigth,length;
	double min_dist;
	int point_count = 30;

	min_dist = 50;
	width = 200;
	heigth = 100;
	length = 50;
	init(width, heigth,length, min_dist);
	process_list(width, heigth, length, min_dist, point_count);
	create_spheres_pbrt();

	return 0;
}