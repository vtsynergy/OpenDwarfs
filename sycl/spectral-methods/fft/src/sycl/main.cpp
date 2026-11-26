// main.cpp -- option handling of opencl/spectral-methods/fft/src/opencl/common/main.cpp (SHOC driver)
#include <iostream>
#include "dwarf_sycl.hpp"
#include "OptionParser.h"
#include "InvalidArgValue.h"

void addBenchmarkSpecOptions(OptionParser &op);
void RunBenchmark(OptionParser &op);
extern int fft_valid;

int main(int argc, char *argv[])
{
	dsycl::parse_args(&argc, &argv);
	int ret = 0;
	try
	{
		OptionParser op;
		op.addOption("platform", OPT_INT, "0", "specify OpenCL platform to use", 'p');
		op.addOption("device", OPT_VECINT, "", "specify device(s) to run on", 'd');
		op.addOption("passes", OPT_INT, "10", "specify number of passes", 'n');
		op.addOption("size", OPT_VECINT, "1", "specify problem size", 's');
		op.addOption("infoDevices", OPT_BOOL, "", "show info for available platforms and devices", 'i');
		op.addOption("verbose", OPT_BOOL, "", "enable verbose output", 'v');
		op.addOption("quiet", OPT_BOOL, "", "write minimum necessary to standard output", 'q');
		addBenchmarkSpecOptions(op);
		if (!op.parse(argc, argv))
		{
			op.usage();
			return (op.HelpRequested() ? 0 : 1 );
		}
		if (op.getOptionBool("infoDevices")) { std::cerr << "use sycl-ls for device information\n"; return 0; }
		RunBenchmark(op);
	}
	catch( InvalidArgValue& e ) { std::cerr << e.what() << ": " << e.GetMessage() << std::endl; ret = 1; }
	catch( std::exception& e ) { std::cerr << e.what() << std::endl; ret = 1; }
	catch( ... ) { std::cerr << "unrecognized exception caught" << std::endl; ret = 1; }
	if (ret == 0) return dwarf_report(fft_valid);
	return ret;
}
