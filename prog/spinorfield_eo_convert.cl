//eoprec operations
__kernel void convert_to_eoprec(__global spinorfield_eoprec* even, __global spinorfield_eoprec* odd, __global spinorfield* in)
{
	int spacepos, timepos;
	for(int n=0; n<VOL4D/2; n++) {
		     get_even_site(n, &spacepos, &timepos);
		     even[n] = in[get_global_pos(spacepos,timepos)];
		     get_odd_site(n, &spacepos, &timepos);
		     odd[n] = in[get_global_pos(spacepos,timepos)];
		}
	return;
}

__kernel void convert_from_eoprec(__global spinorfield_eoprec* even, __global spinorfield_eoprec* odd, __global spinorfield* out){
  int id = get_global_id(0);
  if(id ==0){
	int spacepos, timepos;
	for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
		get_even_site(n, &spacepos, &timepos);
		out[get_global_pos(spacepos,timepos)] = even[n];
		get_odd_site(n, &spacepos, &timepos);
		out[get_global_pos(spacepos,timepos)] = odd[n];
	}
   }
	return;
}
