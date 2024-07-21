#define MPC_DEBUG
#define MPC_IMPLEMENTATION
#include "mpc.h"
#include <stdio.h>

int main() {
	Program program = mpc_load_program("meta.c"); (void)program;
	//printf("%d\n", program.size);
	//printf("%s\n", program.bytes);
	//Token my_type = {"my_type", 8};
	
	//mpc_token_recognize(program, my_type);
	//int n = mpc_write_program("copy.h", program);
	//printf("%d\n", n);
	return 0;
}
