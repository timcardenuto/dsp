#include "generateGenericPRN.h"


static int polynomial_order = 0;
static std::string initial_fill_str;
static std::string feedback_taps_str;
static std::string output_taps_str;
static std::vector<int> initial_fill;
static std::vector<int> feedback_taps;
static std::vector<int> output_taps;
static int verbose = 0;
static int vv = 0;
std::string delimiter = ",";


/* help message block */
void displayCmdUsage() {
	puts("Usage: ./generateGenericPRN [OPTIONS] \n\
	-p	--polynomial-order		Plots the geolocations, target locations and ellipses \n\
	-i	--initial-fill		Uses shared memory kernels for comparison \n\
	-f	--feedback-taps		Uses constant memory kernels for comparison\n\
	-o	--output-taps	Path to configuration file. Any parameters specified on the command line \n\
				will override the equivalent ones in this file \n\
	-v	--verbose	Prints additional output \n\
		--vv		Additional debug type information, including the output of EVERY operation \n\
		--help		Display this message \n");
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
			{"polynomial-order", required_argument, 	0, 'p'},
			{"initial-fill",     required_argument, 	0, 'i'},
			{"feedback-taps",    required_argument, 	0, 'f'},
			{"output-taps",      required_argument, 	0, 'o'},
			{"vv",               no_argument, 	      &vv, 1},
			{"verbose",          no_argument, 	 &verbose, 1},
			{"help",             no_argument, 		    0, 'h'},
			{0,                  0, 					0, 0},
		};

		c = getopt_long_only(argc, argv, "p:i:f:o:hv", long_options, &option_index);

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
			case 'p':
				polynomial_order = strtoul(optarg, &ptr, 10);
				if (strcmp(ptr,"")) {
					printf("Value %s of option %s is not a number \n", ptr, long_options[option_index].name);
					exit(1);
				}
		        break;
			case 'i':
				initial_fill_str = optarg;
				initial_fill = parseVectorString(initial_fill_str);
				break;
       		case 'f':
				feedback_taps_str = optarg;
				feedback_taps = parseVectorString(feedback_taps_str);
		        break;
       		case 'o':
				output_taps_str = optarg;
				output_taps = parseVectorString(output_taps_str);
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
		for (int j=order-1; j>=0; j--) {
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
			if (j != 0) { fill[j] = fill[j-1]; }
		}
		fill[0] = feedback;
	}
	return output;
}


int main(int argc, char **argv) {
    printf("##### Generating Psuedo-Random Number Sequence\n");

    checkCmdArgs(argc, argv);

	// print out args
	std::cout << "Polynomial Order: " << polynomial_order << std::endl;
	std::cout << "Initial Fill:     ";
	for (int i=0; i<initial_fill.size(); i++) {
		std::cout << initial_fill[i] << ", ";
	}
	std::cout << "\nFeedback Taps:    ";
	for (int i=0; i<feedback_taps.size(); i++) {
		std::cout << feedback_taps[i] << ", ";
	}
	std::cout << "\nOutput Taps:      ";
	for (int i=0; i<output_taps.size(); i++) {
		std::cout << output_taps[i] << ", ";
	}
	std::cout << "\n" << std::endl;

	//std::cout << convertVectorToDecimal(initial_fill) << std::endl;

	std::vector<int> PRN = generatePRN(polynomial_order, initial_fill, feedback_taps, output_taps);

	std::cout << "\nPRN:      ";
	for (int i=0; i<PRN.size(); i++) {
		std::cout << PRN[i] << ", ";
	}
	std::cout << "\n" << std::endl;
}
