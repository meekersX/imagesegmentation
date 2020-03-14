#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cfloat>
#include <Magick++.h>
#include "luv.h"
#include "kmeans.h"
using namespace std;
using namespace Magick;

vector<double> featurevector(unsigned int start_x, unsigned int start_y, const Image &image);
pair<vector<vector<double> >, vector<unsigned int> > optimal(const vector<vector<double> > &fvs);

vector<double> premeans(6, 0.0);
vector<double> prestddevs(6, 1.0);
vector<double> maximums(6, 1.0);
vector<double> minimums(6, 1.0);

void printvector(const vector<double> &fv, const string header)
{
	unsigned int  i;

	cout << header << ":";
	for (i = 0; i < fv.size(); i++)
	{
		cout << " " << fv[i];
	}
	cout << endl;
}

void scale(vector<vector<double> > &fvs)
{
	unsigned int i;
	unsigned int j;

	premeans.resize(6);
	prestddevs.resize(6);
	maximums.resize(6);
	minimums.resize(6);
	fill(premeans.begin(), premeans.end(), 0.0);
	fill(prestddevs.begin(), prestddevs.end(), 0.0);
	fill(maximums.begin(), maximums.end(), DBL_MIN);
	fill(minimums.begin(), minimums.end(), DBL_MAX);

	// calculate premeans
	for (i = 0; i < fvs.size(); i++)
		for (j = 0; j < 6; j++)
			premeans[j] += fvs[i][j];

	for (j = 0; j < 6; j++)
		premeans[j] /= fvs.size();

	printvector(premeans, "premeans");

	// make means 0.0
	for (i = 0; i < fvs.size(); i++)
		for (j = 0; j < 6; j++)
			fvs[i][j] -= premeans[j];

	// calculate prestddevs
	for (i = 0; i < fvs.size(); i++)
		for (j = 0; j < 6; j++)
			prestddevs[j] += fvs[i][j]*fvs[i][j];

	for (i = 0; i < 6; i++)
		prestddevs[i] = sqrt(prestddevs[i] / fvs.size());

	printvector(prestddevs, "prestddevs");

	// make stdevs 1.0
	for (i = 0; i < fvs.size(); i++)
		for (j = 0; j < 6; j++)
			fvs[i][j] /= prestddevs[j];

	// calculate maximums and minimums
	for (i = 0; i < fvs.size(); i++)
		for (j = 0; j < 6; j++)
		{
			if (fvs[i][j] > maximums[j])
				maximums[j] = fvs[i][j];
			if (fvs[i][j] < minimums[j])
				minimums[j] = fvs[i][j];
		}

	// scale by maximums on minimums
	for (i = 0; i < fvs.size(); i++)
		for (j = 0; j < 6; j++)
		{
			fvs[i][j] -= (maximums[j] + minimums[j]) / 2;
			fvs[i][j] /= (maximums[j] - minimums[j]);
		}
}

void antiscale(vector<double> &fv)
{
	unsigned int i;

	for (i = 0; i < 6; i++)
	{
		fv[i] *= (maximums[i] - minimums[i]);
		fv[i] += (maximums[i] + minimums[i]) / 2;
		fv[i] *= prestddevs[i];
		fv[i] += premeans[i];
	}
}

bool big(const pair<unsigned int, unsigned int> &a, const pair<unsigned int, unsigned int> &b)
{
	return (a.second > b.second);
}

bool small(const pair<unsigned int, unsigned int> &a, const pair<unsigned int, unsigned int> &b)
{
	return (a.second < b.second);
}

double calculatedistance(const vector<double> &a, const vector<double> &b)
{
	double distance = 0.0;
	unsigned int i;

	for (i = 0; i < 6; i++)
		distance += (a[i] - b[i])*(a[i] - b[i]);

	return distance;
}

double distancetomeans(const vector<double> &a, const vector<vector<double> > &means, unsigned int n)
{
	unsigned int i;
	double distance = 0.0;

	for (i = 0; i < n; i++)
		distance += calculatedistance(a, means[i]);

	return distance;
}

void writeregions(const vector<vector<double> > &fvs, vector<vector<double> > means, const vector<unsigned int> &sets, unsigned int width, unsigned int height)
{
	Image regions(Geometry(width / 4, height / 4), Color(0, 0, 0));
	vector<Color> colors(means.size());
	unsigned int i;
	unsigned int x;
	unsigned int y;
	double R;
	double G;
	double B;

	for (i = 0; i < means.size(); i++)
	{
		antiscale(means[i]);
		LUVtoRGB(means[i][0], means[i][1], means[i][2], R, G, B);
		colors[i] = ColorRGB(R, G, B);
	}

	i = 0;
	for (y = 0; y < height / 4; y++)
	{
		for (x = 0; x < width / 4; x++)
		{
			regions.pixelColor(x, y, colors[sets[i]]);
			i++;
		}
	}

	regions.write("regions.gif");
}

void writefvimages(vector<vector<double> > fvs, unsigned int width, unsigned int height)
{
	Image luv(Geometry(width / 2, height / 2), Color(0, 0, 0));
	Image daub(Geometry(width / 2, height / 2), Color(0, 0, 0));
	vector<double> tempfv(6);
	vector<double> maximum(6, DBL_MIN);
	vector<double> minimum(6, DBL_MAX);
	unsigned int i;
	unsigned int j;
	unsigned int x;
	unsigned int y;

	for (i = 0; i < fvs.size(); i++)
		antiscale(fvs[i]);

	i = 0;
	for (y = 0; y < height / 4; y++)
	{
		for (x = 0; x < width / 4; x++)
		{
			tempfv = fvs[i];
			LUVtoRGB(tempfv[0], tempfv[1], tempfv[2], tempfv[0], tempfv[1], tempfv[2]);
			ColorRGB color(tempfv[0], tempfv[1], tempfv[2]);
			luv.pixelColor(x, y, color);
			daub.pixelColor(x, y, color);
			for (j = 0; j < 6; j++)
			{
				if (fvs[i][j] > maximum[j])
					maximum[j] = fvs[i][j];
				if (fvs[i][j] < minimum[j])
					minimum[j] = fvs[i][j];
			}
			i++;
		}
	}

	for (i = 0; i < fvs.size(); i++)
		for (j = 0; j < 6; j++)
		{
			fvs[i][j] -= minimum[j];
			fvs[i][j] /= (maximum[j] - minimum[j]);
		}

	i = 0;
	for (y = 0; y < height / 4; y++)
	{
		for (x = width / 4; x < width / 2; x++)
		{
			luv.pixelColor(x, y, ColorRGB(fvs[i][0], fvs[i][0], fvs[i][0]));
			daub.pixelColor(x, y, ColorRGB(fvs[i][3], fvs[i][3], fvs[i][3]));
			i++;
		}
	}

	i = 0;
	for (y = height / 4; y < height / 2; y++)
	{
		for (x = 0; x < width / 4; x++)
		{
			luv.pixelColor(x, y, ColorRGB(fvs[i][1], fvs[i][1], fvs[i][1]));
			daub.pixelColor(x, y, ColorRGB(fvs[i][4], fvs[i][4], fvs[i][4]));
			i++;
		}
	}

	i = 0;
	for (y = height / 4; y < height / 2; y++)
	{
		for (x = width / 4; x < width / 2; x++)
		{
			luv.pixelColor(x, y, ColorRGB(fvs[i][2], fvs[i][2], fvs[i][2]));
			daub.pixelColor(x, y, ColorRGB(fvs[i][5], fvs[i][5], fvs[i][5]));
			i++;
		}
	}

	luv.write("luv.png");
	daub.write("daub.png");
}

unsigned int best(const vector<unsigned int> &sets, unsigned int width, unsigned int height, unsigned int x, unsigned int y)
{
	vector<unsigned int> count_sets;
	unsigned int best;
	unsigned int max;
	unsigned int current;
	unsigned int i;

	count_sets.reserve(8);

	if (x - 1 >= 0)
		count_sets.push_back(sets[y*width + x - 1]);

	if (x + 1 < width)
		count_sets.push_back(sets[y*width + x + 1]);

	if (y - 1 >= 0)
		count_sets.push_back(sets[(y - 1)*width + x]);

	if (y + 1 < height)
		count_sets.push_back(sets[(y + 1)*width + x]);

	if (x - 1 >= 0 && y - 1 >= 0)
		count_sets.push_back(sets[(y - 1)*width + x - 1]);

	if (x - 1 >= 0 && y + 1 >= 0)
		count_sets.push_back(sets[(y + 1)*width + x - 1]);

	if (x + 1 >= 0 && y - 1 >= 0)
		count_sets.push_back(sets[(y - 1)*width + x + 1]);

	if (x + 1 >= 0 && y + 1 >= 0)
		count_sets.push_back(sets[(y + 1)*width + x + 1]);

	sort(count_sets.begin(), count_sets.end());

	for (i = 0; count_sets[i] == count_sets[0]; i++);

	best = count_sets[0];
	max = i;

	current = 0;
	for (i++; i < count_sets.size(); i++)
	{
		current++;
		if (count_sets[i] != count_sets[i - 1])
		{
			if (current > max)
			{
				max = current;
				best = count_sets[i - 1];
			}
			current = 0;
		}
	}

	if (current > 0)
	{
		if (current > max)
		{
			max = current;
			best = count_sets[i - 1];
		}
	}

	return best;
}

bool and_set(const vector<unsigned int> &sets, unsigned int width, unsigned int height, signed int x, signed int y)
{
	bool result = true;
	unsigned int set = sets[y*width + x];

	if (x - 1 >= 0)
		result = result && set == sets[y*width + x - 1];

	if (x + 1 < width)
		result = result && set == sets[y*width + x + 1];

	if (y - 1 >= 0)
		result = result && set == sets[(y - 1)*width + x];

	if (y + 1 < height)
		result = result && set == sets[(y + 1)*width + x];
/*
	if (x - 1 >= 0 && y - 1 >= 0)
		result = result && set == sets[(y - 1)*width + x - 1];

	if (x - 1 >= 0 && y + 1 >= 0)
		result = result && set == sets[(y + 1)*width + x - 1];

	if (x + 1 >= 0 && y - 1 >= 0)
		result = result && set == sets[(y - 1)*width + x + 1];

	if (x + 1 >= 0 && y + 1 >= 0)
		result = result && set == sets[(y + 1)*width + x + 1];
*/
	return result;
}

bool or_set(const vector<unsigned int> &sets, unsigned int width, unsigned int height, signed int x, signed int y)
{
	bool result = false;
	unsigned int set = sets[y*width + x];

	if (x - 1 >= 0)
		result = result || set == sets[y*width + x - 1];

	if (x + 1 < width)
		result = result || set == sets[y*width + x + 1];

	if (y - 1 >= 0)
		result = result || set == sets[(y - 1)*width + x];

	if (y + 1 < height)
		result = result || set == sets[(y + 1)*width + x];
/*
	if (x - 1 >= 0 && y - 1 >= 0)
		result = result || set == sets[(y - 1)*width + x - 1];

	if (x - 1 >= 0 && y + 1 >= 0)
		result = result || set == sets[(y + 1)*width + x - 1];

	if (x + 1 >= 0 && y - 1 >= 0)
		result = result || set == sets[(y - 1)*width + x + 1];

	if (x + 1 >= 0 && y + 1 >= 0)
		result = result || set == sets[(y + 1)*width + x + 1];
*/
	return result;
}

vector<unsigned int> erode(const vector<unsigned int> &sets, unsigned int width, unsigned int height)
{
	unsigned int x;
	unsigned int y;
	bool result;
	vector<unsigned int> new_sets(sets.size());

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++)
		{
			result = and_set(sets, width, height, x, y);

			if (result)
			{
				new_sets[y*width + x] = sets[y*width + x];
			}
			else
			{
				new_sets[y*width + x] = best(sets, width, height, x, y);
			}
		}

	return new_sets;
}

vector<unsigned int> dilate(const vector<unsigned int> &sets, unsigned int width, unsigned int height)
{
	unsigned int x;
	unsigned int y;
	bool result;
	vector<unsigned int> new_sets(sets.size());

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++)
		{
			result = or_set(sets, width, height, x, y);

			if (result)
			{
				new_sets[y*width + x] = sets[y*width + x];
			}
			else
			{
				new_sets[y*width + x] = best(sets, width, height, x, y);
			}
		}

	return new_sets;
}

int main(int argc, char** argv)
{
	InitializeMagick(*argv);
    
	Image image("working.jpg");
	Image correct("211.gif");
	vector<vector<double> > fvs;
	vector<unsigned int> sets;
	vector<unsigned int> correct_sets;
	vector<vector<double> > means;
	vector<unsigned int> sizes;
	pair<vector<vector<double> >, vector<unsigned int> > value;
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int x;
	unsigned int y;
	map<long long, unsigned int> set_mapping;
	long long raw;
	vector<pair<unsigned int, unsigned int> > correct_mapping;
	vector<pair<unsigned int, unsigned int> > current_mapping;
	vector<unsigned int> current_map;
	vector<vector<double> > population(1000, vector<double>(4));
	vector<vector<double> > fit(population.size() * 0.1, vector<double>(4));
	vector<pair<unsigned int, unsigned int> > population_fitness(population.size());
	unsigned int different;
	double sum;

//	srand(2);
	srand(time(NULL));
/*
	double R, G, B, X, Y, Z, L, U, V;
	for (R = 0.0; R <= 1.0; R += 1.0)
		for (G = 0.0; G <= 1.0; G += 1.0)
			for (B = 0.0; B <= 1.0; B += 1.0)
			{
				cout << "R: " << R << " G: " << G << " B: " << B << endl;
				RGBtoXYZ(R, G, B, X, Y, Z);
				cout << "X: " << X << " Y: " << Y << " Z: " << Z << endl;
				XYZtoLUV(X, Y, Z, L, U, V);
				cout << "L: " << L << " U: " << U << " V: " << V << endl;
				cout << endl;
			}
	return 0;
*/
	// calculate set information for correct segmentation
	for (y = 0; y < correct.rows(); y++)
	{
		for (x = 0; x < correct.columns(); x++)
		{
			ColorRGB pixel = correct.pixelColor(x, y);
			raw = ((int)(8*pixel.red()) << 16) | ((int)(8*pixel.green()) << 8) | ((int)(8*pixel.blue()));

			if (set_mapping.find(raw) == set_mapping.end())
			{
				set_mapping[raw] = set_mapping.size();
				correct_mapping.push_back(pair<unsigned int, unsigned int>(set_mapping.size() - 1, 0));
			}

			correct_sets.push_back(set_mapping[raw]);
			correct_mapping[set_mapping[raw]].second++;
		}
	}
	sort(correct_mapping.begin(), correct_mapping.end(), big);

//	for (i = 0; i < correct_mapping.size(); i++)
//		cout << correct_mapping[i].first << ": " << correct_mapping[i].second << endl;

	// calculate fv's for image
	for (y = 0; y < image.rows(); y += 4)
	{
		for (x = 0; x < image.columns(); x += 4)
		{
			fvs.push_back(featurevector(x, y, image));
		}
	}

	scale(fvs);
//	printvector(calcscales(fvs), "scales");

	writefvimages(fvs, image.columns(), image.rows());

	value = optimal(fvs);
	means = value.first;
	sets = value.second;
/*
	fvs.resize(4*fvs.size());
	sets.resize(4*sets.size());

	for (y = image.rows() / 2 - 1; (signed)y >= 0; y--)
		for (x = image.columns() / 2 - 1; (signed)x >= 0; x--)
		{
			fvs[y*image.columns()/2+x] = fvs[y/2*image.columns()/4+x/2];
			sets[y*image.columns()/2+x] = sets[y/2*image.columns()/4+x/2];
//			cout << y*image.columns()/2+x << " -> " << y/2*image.columns()/4+x/2 << endl;
		}
*/

	sets = dilate(sets, image.columns() / 4, image.rows() / 4);
	sets = erode(sets, image.columns() / 4, image.rows() / 4);
	sets = erode(sets, image.columns() / 4, image.rows() / 4);
	sets = dilate(sets, image.columns() / 4, image.rows() / 4);

	sizes.resize(means.size());

	for (i = 0; i < means.size(); i++)
	{
		sizes[i] = 0;
		for (j = 0; j < 6; j++)
		{
			means[i][j] = 0.0;
		}
	}

	for (i = 0; i < fvs.size(); i++)
	{
		int test = sets[i];
		sizes[sets[i]]++;
		for (j = 0; j < 6; j++)
		{
			means[sets[i]][j] += fvs[i][j];
		}
	}

	for (i = 0; i < means.size(); i++)
	{
		for (j = 0; j < 6; j++)
		{
			means[i][j] /= sizes[i];
		}
	}

	writeregions(fvs, means, sets, image.columns(), image.rows());

/*
	for (i = 0; i < population.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			population[i][j] = rand() / (double)RAND_MAX;
		}
	}

	while (true)
	{
		for (i = 0; i < population.size(); i++)
		{
			weights = population[i];

			value = optimal(fvs);
			means = value.first;
			sets = value.second;

			current_mapping.resize(means.size());
			for (j = 0; j < means.size(); j++)
			{
				current_mapping[j].first = j;
				current_mapping[j].second = 0;
			}

			for (j = 0; j < sets.size(); j++)
				current_mapping[sets[j]].second++;

			sort(current_mapping.begin(), current_mapping.end(), big);

	//		for (j = 0; j < current_mapping.size(); j++)
	//			cout << current_mapping[j].first << ": " << current_mapping[j].second << endl;

			current_map.resize(means.size());
			for (j = 0; j < means.size(); j++)
				current_map[current_mapping[j].first] = correct_mapping[j].first;

	//		for (j = 0; j < current_map.size(); j++)
	//			cout << j << "->" << current_map[j] << endl;

			different = 0;
			for (j = 0; j < sets.size(); j++)
				if (current_map[sets[j]] != correct_sets[j])
					different++;

	//		cout << "Different: " << different << endl;
			cout << ".";
			if ((i + 1) % 50 == 0)
				cout << endl;

			population_fitness[i].first = i;
			population_fitness[i].second = different;
		}

		sum = 0;
		for (i = 0; i < population_fitness.size(); i++)
			sum += population_fitness[i].second;

		cout << "Average: " << sum / population_fitness.size();

		sort(population_fitness.begin(), population_fitness.end(), small);

		sum = 0;
		for (i = 0; i < fit.size(); i++)
			sum += population_fitness[i].second;

		cout << "\tTop Average: " << sum / fit.size() << endl;

		// from top 10% of population make new population by breeding
		for (i = 0; i < fit.size(); i++)
		{
			for (j = 0; j < 4; j++)
			{
				fit[i][j] = population[population_fitness[i].first][j];
			}			cout << (i + 1) << " (" << population_fitness[i].second << ")";
			printvector(fit[i], "");
		}

		x = 0;
		do
		{
			for (i = 0; i < fit.size(); i++)
				for (j = i + 1; j < fit.size() && x < population.size(); j++, x++)
					for (k = 0; k < 4; k++)
					{
						switch (rand() % 2)
						{
						case 0:
							population[x][k] = fit[i][k];
							break;
						case 1:
							population[x][k] = fit[j][k];
							break;
						}
//						population[x][k] = population[rand() % 10][k];
//						cout << x << "/" << k << "=" << population[x][k] << endl;
					}
//		printvector(population[0], "Breeded");
		} while (x < population.size());

		// mutate new population
		for (i = 0; i < population.size(); i++)
		{
//			unsigned int selected = rand() % population.size();
			unsigned int selected = i;
			unsigned int feature = rand() % 4;

//			cout << "Old: " << population[selected][feature];

//			if (rand() % 2)
//				population[selected][feature] *= 1.0 - 2.0 * log(rand() / (double)RAND_MAX);
//			else
//				population[selected][feature] *= rand() / (double)RAND_MAX;

			population[selected][feature] *= rand() / (double)RAND_MAX * 2.0;
//			cout << "\tNew: " << population[selected][feature] << endl;
		}

		weights = fit[0];

		value = optimal(fvs);
		means = value.first;
		sets = value.second;

		writeregions(fvs, means, sets, image.columns(), image.rows());
	}
*/
	current_mapping.resize(means.size());
	for (j = 0; j < means.size(); j++)
	{
		current_mapping[j].first = j;
		current_mapping[j].second = 0;
	}

	for (j = 0; j < sets.size(); j++)
		current_mapping[sets[j]].second++;

	sort(current_mapping.begin(), current_mapping.end(), big);

	cout << "Correct: " << endl;
	for (j = 0; j < correct_mapping.size(); j++)
		cout << correct_mapping[j].first << ": " << correct_mapping[j].second << endl;

	cout << "Current: " << endl;
	for (j = 0; j < current_mapping.size(); j++)
		cout << current_mapping[j].first << ": " << current_mapping[j].second << endl;

	current_map.resize(means.size());
	for (j = 0; j < means.size(); j++)
		current_map[current_mapping[j].first] = correct_mapping[j].first;

	for (j = 0; j < current_map.size(); j++)
		cout << j << "->" << current_map[j] << endl;

	different = 0;
	for (j = 0; j < sets.size(); j++)
	{
		if (current_map[sets[j]] != correct_sets[j])
		{
			different++;
		}
//		cout << sets[j] << " " << current_map[sets[j]] << " " << correct_sets[j] << endl;
//		break;
	}
	cout << "Different: " << different << endl;

	return 0;
}

vector<double> featurevector(unsigned int start_x, unsigned int start_y, const Image &image)
{
	vector<double> fv(6);
	ColorRGB color;
	unsigned int x;
	unsigned int y;
	unsigned int i;
	double t[4][4];
	double L;
	double U;
	double V;

	fv[0] = 0.0;
	fv[1] = 0.0;
	fv[2] = 0.0;

	for (x = 0; x < 4; x++)
	{
		for (y = 0; y < 4; y++)
		{
			color = image.pixelColor(start_x + x, start_y + y);
			RGBtoLUV(color.red(), color.green(), color.blue(), L, U, V);
			fv[0] += L;
			fv[1] += U;
			fv[2] += V;
			t[y][x] = L;
		}
	}

	// Y
	fv[0] /= 16;
	// U
	fv[1] /= 16;
	// V
	fv[2] /= 16;

	// Daubechies-4 Transform (lifting)
	// s = scaling = lowpass  = L
	// d = wavelet = highpass = H

	// t[i][0] =  L[i][0]
	// t[i][1] = -H[i][1]
	// t[i][2] =  L[i][1]
	// t[i][3] = -H[i][0]

	for (i = 0; i < 4; i++)
	{
		t[i][0] += sqrt(3.0) * t[i][1];
		t[i][2] += sqrt(3.0) * t[i][3];

		t[i][1] -= sqrt(3.0)/4 * t[i][0] + (sqrt(3.0)-2)/4 * t[i][2];
		t[i][3] -= sqrt(3.0)/4 * t[i][2] + (sqrt(3.0)-2)/4 * t[i][0];

		t[i][0] -= t[i][3];
		t[i][2] -= t[i][1];

		t[i][0] *= (sqrt(3.0)-1) / sqrt(2.0);
		t[i][1] *= (sqrt(3.0)+1) / sqrt(2.0);
		t[i][2] *= (sqrt(3.0)-1) / sqrt(2.0);
		t[i][3] *= (sqrt(3.0)+1) / sqrt(2.0);
	}

	for (i = 0; i < 4; i++)
	{
		t[0][i] += sqrt(3.0) * t[1][i];
		t[2][i] += sqrt(3.0) * t[3][i];

		t[1][i] -= sqrt(3.0)/4 * t[0][i] + (sqrt(3.0)-2)/4 * t[2][i];
		t[3][i] -= sqrt(3.0)/4 * t[2][i] + (sqrt(3.0)-2)/4 * t[0][i];

		t[0][i] -= t[3][i];
		t[2][i] -= t[1][i];

		t[0][i] *= (sqrt(3.0)-1) / sqrt(2.0);
		t[1][i] *= (sqrt(3.0)+1) / sqrt(2.0);
		t[2][i] *= (sqrt(3.0)-1) / sqrt(2.0);
		t[3][i] *= (sqrt(3.0)+1) / sqrt(2.0);
	}

	// HL
	fv[3] = sqrt(t[0][1]*t[0][1] + t[0][3]*t[0][3] + t[2][1]*t[2][1] + t[2][3]*t[2][3]) / 2;
	// LH
	fv[4] = sqrt(t[1][0]*t[1][0] + t[1][2]*t[1][2] + t[3][0]*t[3][0] + t[3][2]*t[3][2]) / 2;
	// HH
	fv[5] = sqrt(t[1][1]*t[1][1] + t[1][3]*t[1][3] + t[3][1]*t[3][1] + t[3][3]*t[3][3]) / 2;

/*
	ofstream out("debug.txt", ofstream::out | ofstream::app);
	out << fv[0];
	out << "\t" << fv[1];
	out << "\t" << fv[2];
	out << "\t" << fv[3];
	out << "\t" << fv[4];
	out << "\t" << fv[5];
	out << endl;
*/
	return fv;
}

pair<vector<vector<double> >, vector<unsigned int> > optimal(const vector<vector<double> > &fvs)
{
	unsigned int k = 6;
	vector<vector<double> > means(k, vector<double>(6, 0.0));
	double distance;
	double max_distance;
	double D;
	unsigned int i;
	unsigned int j;
	unsigned int farthest;

	for (i = 0; i < fvs.size(); i++)
	{
		for (j = 0; j < 6; j++)
		{
			means[0][j] += fvs[i][j];
		}
	}

//	cout << endl << "Mean 1:";
	for (j = 0; j < 6; j++)
	{
		means[0][j] /= fvs.size();
//		cout << " " << means[0][j];
	}
//	cout << endl;

	for (i = 1; i < k; i++)
	{
		max_distance = distancetomeans(fvs[0], means, i);
		farthest = 0;

		for (j = 1; j < fvs.size(); j++)
		{
			distance = distancetomeans(fvs[j], means, i);

			if (distance > max_distance)
			{
				max_distance = distance;
				farthest = j;
			}
		}

//		cout << "Mean " << (i + 1) << " (" << farthest << "):";
		for (j = 0; j < 6; j++)
		{
			means[i][j] = fvs[farthest][j];
//			cout << " " << means[i][j];
		}
//		cout << endl;
	}

	return kmeans(fvs, means, k, D);
}