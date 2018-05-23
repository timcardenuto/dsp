#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>		// for getopt support
#include <getopt.h>		// for getopt_long support
#include <math.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector


void displayCmdUsage();
std::vector<int> parseVectorString(std::string s);
int convertVectorToDecimal(std::vector<int> v);
void checkCmdArgs(int argc, char **argv);
std::vector<int> generatePRN(int order, std::vector<int> fill, std::vector<int> ftaps, std::vector<int> otaps);
int main(int argc, char **argv);
