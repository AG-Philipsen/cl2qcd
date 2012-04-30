//If the include does not work then the contents of the ranluxcl.cl file can just be
//pasted here at the top instead.
#include "ranluxcl.cl"

__kernel void Kernel_Ranluxcl_Init(
	private int ins,
	global float4 *ranluxcltab)
{
	ranluxcl_initialization(ins, ranluxcltab);
}

__kernel void Kernel_PRN(private int KernelCycles,
                         global float4* ranluxcltab,
                         global float* PRNs)
{
	//Downloading ranluxcltab. The state of RANLUXCL is stored in ranluxclstate.
	ranluxcl_state_t ranluxclstate;
	ranluxcl_download_seed(&ranluxclstate, ranluxcltab);

	float4 randomnr;

	//Generate some numbers
	for(int i=0; i<KernelCycles; i+=4)
		randomnr = ranluxcl(&ranluxclstate);

	//Uploading only last number generated.
	PRNs[get_global_id(0)] = randomnr.w;

	//Uploading ranluxcltab
	ranluxcl_upload_seed(&ranluxclstate, ranluxcltab);
}
