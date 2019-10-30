#include "picojson.h"
#include "metadata.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;
using namespace picojson;


struct metadata md;

int json_parse(const char *inputfile){

	// Declaration
	stringstream ss;
	ifstream	f;
	unsigned int i;

	// Read Json file
	f.open(inputfile,ios::binary);
	if (!f.is_open()) return 1;
	ss << f.rdbuf();
	f.close();

	// Parse Json data
	value v;
	ss >> v;

	string err = get_last_error();
	if (!err.empty()) {
		cerr << err << endl;
		return -1;
	}

	object& image = v.get<object>()["image"].get<object>();
	object& data = v.get<object>()["devices"].get<object>();
	object& mla = data["mla"].get<object>();
	object& sensorOffset = mla["sensorOffset"].get<object>();
	object& sensor = data["sensor"].get<object>();

	md.w = image["width"].get<double>();
	md.h = image["height"].get<double>();
	md.pP = sensor["pixelPitch"].get<double>();
	md.x_off = sensorOffset["x"].get<double>();
	md.y_off = sensorOffset["y"].get<double>();
	md.z_off = sensorOffset["z"].get<double>();
	md.lP = mla["lensPitch"].get<double>();
	md.rot = mla["rotation"].get<double>();
}

/*
int main(int argc,char *argv[])
{
	json_parse(md,argv[1]);
	cout << md.h << "," << md.w <<endl;
	return 0;
}
*/