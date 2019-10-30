#include "picojson.h"
#include "metadata.h"
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>

//metadata md;
using namespace std;

extern metadata md;

int main(int argc,char *argv[]){
    json_parse(argv[1]);
	cout << md.w<< endl;
    return 0;
}