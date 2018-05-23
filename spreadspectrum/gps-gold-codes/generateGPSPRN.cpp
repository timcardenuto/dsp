/*
    g++ -std=c++11 -o generateGPSPRN generateGPSPRN.cpp -I.
*/
#include "generateGPSPRN.h"

std::map<int, std::vector<int>> gpssig;
static int space_vehicle = 0;
static int prn_num = 0;
static std::string taps_str;
static std::vector<int> G2_output_taps = {0, 0};
static int verbose = 0;
static int vv = 0;
std::string delimiter = ",";


void buildSVMap() {
    gpssig.insert(std::make_pair(1, std::vector<int>{2,6}));
    gpssig.insert(std::make_pair(2, std::vector<int>{3,7}));
    gpssig.insert(std::make_pair(3, std::vector<int>{4,8}));
    gpssig.insert(std::make_pair(4, std::vector<int>{5,9}));
    gpssig.insert(std::make_pair(5, std::vector<int>{1,9}));
    gpssig.insert(std::make_pair(6, std::vector<int>{2,10}));
    gpssig.insert(std::make_pair(7, std::vector<int>{1,8}));
    gpssig.insert(std::make_pair(8, std::vector<int>{2,9}));
    gpssig.insert(std::make_pair(9, std::vector<int>{3,10}));
    gpssig.insert(std::make_pair(10, std::vector<int>{2,3}));
    gpssig.insert(std::make_pair(11, std::vector<int>{3,4}));
    gpssig.insert(std::make_pair(12, std::vector<int>{5,6}));
    gpssig.insert(std::make_pair(13, std::vector<int>{6,7}));
    gpssig.insert(std::make_pair(14, std::vector<int>{7,8}));
    gpssig.insert(std::make_pair(15, std::vector<int>{8,9}));
    gpssig.insert(std::make_pair(16, std::vector<int>{9,10}));
    gpssig.insert(std::make_pair(17, std::vector<int>{1,4}));
    gpssig.insert(std::make_pair(18, std::vector<int>{2,5}));
    gpssig.insert(std::make_pair(19, std::vector<int>{3,6}));
    gpssig.insert(std::make_pair(20, std::vector<int>{4,7}));
    gpssig.insert(std::make_pair(21, std::vector<int>{5,8}));
    gpssig.insert(std::make_pair(22, std::vector<int>{6,9}));
    gpssig.insert(std::make_pair(23, std::vector<int>{1,3}));
    gpssig.insert(std::make_pair(24, std::vector<int>{4,6}));
    gpssig.insert(std::make_pair(25, std::vector<int>{5,7}));
    gpssig.insert(std::make_pair(26, std::vector<int>{6,8}));
    gpssig.insert(std::make_pair(27, std::vector<int>{7,9}));
    gpssig.insert(std::make_pair(28, std::vector<int>{8,10}));
    gpssig.insert(std::make_pair(29, std::vector<int>{1,6}));
    gpssig.insert(std::make_pair(30, std::vector<int>{2,7}));
    gpssig.insert(std::make_pair(31, std::vector<int>{3,8}));
    gpssig.insert(std::make_pair(32, std::vector<int>{4,9}));
}

/* help message block */
void displayCmdUsage() {
	puts("Usage: ./generateGPSPRN [OPTIONS] \n\
	-s	--space-vehicle    Plots the geolocations, target locations and ellipses \n\
	-p	--prn              Uses shared memory kernels for comparison \n\
	-t	--taps             Uses constant memory kernels for comparison\n\
	-v	--verbose          Prints additional output \n\
		--vv               Additional debug type information, including the output of EVERY operation \n\
		--help             Display this message \n");
    std::cout << "\nValid input:" << std::endl;
    buildSVMap();
    for (int i=1; i<=gpssig.size(); i++) {
        std::cout << "        SV: " << i << " with taps: " << gpssig[i][0] << "," << gpssig[i][1] << std::endl;
    }
	exit(1);
}

std::vector<int> parseVectorString(std::string s) {
	size_t pos = 0;
	std::string token;
	std::vector<int> v;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		v.push_back(std::stoi(token));
		s.erase(0, pos + delimiter.length());
	}
	v.push_back(std::stoi(s));
	// sort vector
	std::sort(v.begin(), v.end());
	return v;
}

// MSB -> LSB
int convertVectorToDecimal(std::vector<int> v) {
	int d = 0;
	int j = 0;
	for (int i=v.size()-1; i>=0; i--) {
		if (v[i] == 1) {
			d = d + pow(2.0,j);
		}
		j = j + 1;
	}
	return d;
}

/* Use getopt to read cmd line arguments */
void checkCmdArgs(int argc, char **argv) {
    int c;
	char *ptr;
   	while (1) {
        int option_index = 0;

        static struct option long_options[] = {
			{"space-vehicle", required_argument, 	     0, 's'},
			{"prn",           required_argument, 	     0, 'p'},
			{"taps",          required_argument,         0, 't'},
			{"vv",                  no_argument,       &vv,   1},
			{"verbose",             no_argument,  &verbose,   1},
			{"help",                no_argument,         0, 'h'},
			{0,                               0,         0,   0},
		};

		c = getopt_long_only(argc, argv, "s:p:t:hv", long_options, &option_index);

		if (c == -1) { break; }

       	switch (c) {
			 case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0) {
					break;
				}
				printf ("option %s", long_options[option_index].name);
				if (optarg) {
					printf (" with arg %s", optarg);
				}
				printf ("\n");
				break;
			case 's':
				space_vehicle = strtoul(optarg, &ptr, 10);
				if (strcmp(ptr,"")) {
					printf("Value %s of option %s is not a number \n", ptr, long_options[option_index].name);
					exit(1);
				}
		        break;
			case 'p':
                prn_num = strtoul(optarg, &ptr, 10);
                if (strcmp(ptr,"")) {
                    printf("Value %s of option %s is not a number \n", ptr, long_options[option_index].name);
                    exit(1);
                }
                break;
       		case 't':
				taps_str = optarg;
				G2_output_taps = parseVectorString(taps_str);
		        break;
       		case 'v':
				verbose = 1;
		        break;
       		default:
				displayCmdUsage();
	    }
	}
	if (optind < argc) {
      	printf ("Unrecognized options: ");
      	while (optind < argc) {
        	printf ("%s ", argv[optind++]);
      	}
		printf("\n"); //putchar ('\n');
		displayCmdUsage();
    }
	return;
}

// TODO should first output bit be before or after first shift operation?
std::vector<int> generatePRN(int order, std::vector<int> fill, std::vector<int> ftaps, std::vector<int> otaps) {
	// sequence bit length, also the number of possible orthogonal PRN at this order
	int N = pow(2.0,order) - 1;
	// initialize output sequence with zeros
	std::vector<int> output(N);

	// clock registers N times to produce full repeating sequence
	for (int i=0; i<N; i++) {
		int ftap_index = ftaps.size()-1;
		int otap_index = otaps.size()-1;
		int feedback = 0;
		// loop through all registers from right to left (1 <- highest poly) to check for taps and shift
		for (int j=order-1; j>0; j--) {
            // check if this register is tapped for feedback
            // NOTE don't forget to subtract 1 from tap number, since registers are indexed from 0 not 1
            //      ex. for taps 1 + x3 + x10, or [3, 10], we really look for j = [2, 9]
            if (j == (ftaps[ftap_index]-1)) {
                feedback = feedback ^ fill[j];
                if (ftap_index > 0) { ftap_index--; }
            }
            // check if this register is tapped for output
            if (j == (otaps[otap_index]-1)) {
                output[i] = output[i] ^ fill[j];
                if (otap_index > 0) { otap_index--; }
            }
			// shift bits
			fill[j] = fill[j-1];
		}
		fill[0] = feedback;
	}
	return output;
}


int main(int argc, char **argv) {
    checkCmdArgs(argc, argv);

    if (verbose) { printf("##### Generating GPS Psuedo-Random Number Sequence\n"); }

    buildSVMap();

    if (0 < space_vehicle && space_vehicle < 33) {
        prn_num = space_vehicle;
        G2_output_taps = {gpssig[space_vehicle][0], gpssig[space_vehicle][1]};
    } else if (0 < prn_num && prn_num < 33) {
        space_vehicle = prn_num;
        G2_output_taps = {gpssig[prn_num][0], gpssig[prn_num][1]};
    } else if (0 < G2_output_taps[0] && G2_output_taps[0] < 10) {
        for (int i=1; i <= gpssig.size(); i++) {
            if (gpssig[i][0] == G2_output_taps[0] && gpssig[i][1] == G2_output_taps[1]) {
                space_vehicle = i;
                prn_num = i;
            }
        }
    } else {
        std::cout << "ERROR: No SV, PRN, or G2 taps specified" << std::endl;
        exit(1);
    }

	// print out args
	if (verbose) { std::cout << "\nSpace Vehicle: " << space_vehicle << std::endl; }
    if (verbose) { std::cout << "PRN number:    " << prn_num << std::endl; }
	if (verbose) { std::cout << "G2 taps:       " << G2_output_taps[0] << "," << G2_output_taps[1] << std::endl; }

    // standard inputs for GPS C/A codes
    int polynomial_order = 10;
    std::vector<int> initial_fill(polynomial_order,1);
    std::vector<int> G1_feedback_taps = {3, 10};
    std::vector<int> G1_output_taps = {10};
    std::vector<int> G2_feedback_taps = {2,3,6,8,9,10};

    // generate G1 and G2, G2 is based on user chosen SV, PRN number, or expected taps
    std::vector<int> G1 = generatePRN(polynomial_order, initial_fill, G1_feedback_taps, G1_output_taps);
	std::vector<int> G2 = generatePRN(polynomial_order, initial_fill, G2_feedback_taps, G2_output_taps);

    // xor G1 and G2 for final GPS PRN
    if (verbose) { std::cout << "PRN:" << std::endl; }
    std::cout << "[";
    std::vector<int> PRN(G1.size());
    for (int i=0; i<G1.size()-1; i++) {
        PRN[i] = G1[i] ^ G2[i];
        std::cout << PRN[i] << ",";
    }
    PRN[G1.size()-1] = G1[G1.size()-1] ^ G2[G1.size()-1];
    std::cout << PRN[G1.size()-1] << "]" << std::endl;
}
