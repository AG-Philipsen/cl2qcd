//The geometry-functions assume the ordering
//	x,y,z,t
//of coordinates
//and the assignment
//	mu = 0 -> time
//	mu = 1 -> x
//	mu = 2 -> y
//	mu = 3 -> z
//One can test the self-consistendy here
//However, some tests dont make sense since one would do exactly the same operations twice,
//	hence, getting possibly the same error twice wihtout realizing it.

__kernel void geometry_test()
{

	int id = get_global_id(0);
	if(id > 0) return;
	int x, y, z, t;
	for(x = 0; x < NSPACE; x++) {
		for(y = 0; y < NSPACE; y++) {
			for(z = 0; z < NSPACE; y++) {
				for(t = 0; t < NTIME; t ++) {
					int glob_pos;
					int ns;
					int ns_up;
					int ns_dn;
					int mu;
					int dir;
					uint3 coord;

					coord.x = x;
					coord.y = y;
					coord.z = z;

					ns = get_nspace(coord);
					glob_pos = get_global_pos(ns, t);

					//test get_spacecoord
					uint3 coord_test = get_coord_spatial(glob_pos);
					if(coord.x != coord_test.x || coord.y != coord_test.y || coord.z != coord_test.z) {
						printf("ERROR at get_spacecoord at x = %i y = %i z = %i t = %i\n", x, y, z, t);
						return;
					}

					//test get_allspacecoord
					coord_test = get_allspacecoord(ns);
					if(coord.x != coord_test.x || coord.y != coord_test.y || coord.z != coord_test.z) {
						printf("ERROR at get_allspacecoord at x = %i y = %i z = %i t = %i\n", x, y, z, t);
						return;
					}

					//test get_nspace
					int nspace_test = x + NSPACE * y + NSPACE * NSPACE * z;
					if(nspace_test != get_nspace(coord)) {
						printf("ERROR at get_nspace at x = %i y = %i z = %i t = %i\n", x, y, z, t);
						return;
					}

					//test get_global_pos
					int tmp = ns + VOLSPACE * t;
					if(get_global_pos(ns, t) != tmp) {
						printf("ERROR at get_global_pos at x = %i y = %i z = %i t = %i\n", x, y, z, t);
						return;
					}

					//test get_even/odd_site
					st_index idx;
					if(glob_pos % 2 == 0) {
						idx = get_even_site(glob_pos);
						if(idx.space != ns || idx.time != t) {
							printf("ERROR at get_even_site at x = %i y = %i z = %i t = %i\n", x, y, z, t);
							return;
						}
					} else if(glob_pos % 2 == 1) {
						idx = get_odd_site(glob_pos);
						if(idx.space != ns || idx.time != t) {
							printf("ERROR at get_odd_site at x = %i y = %i z = %i t = %i\n", x, y, z, t);
							return;
						}
					}

					//test eoprec-indices
					//No test implemented yet

					//test neighbour functions for all directions
					//mu = 0
					// for time-direction there are no neighbour-functions (yet)
					mu = 0;

					//mu = 1
					mu = 1;
					dir = 1;
					ns_up = get_neighbor(ns, mu);
					coord.x = (x + 1) % NSPACE;
					if( get_global_pos(ns_up, t) != get_global_pos(get_nspace(coord), t) ) {
						printf("ERROR at mu = %i in direction %i\n", mu, dir);
						return;
					}
					dir = -1;
					ns_dn = get_lower_neighbor(ns, mu);
					coord.x = (x - 1 + NSPACE) % NSPACE;
					if( get_global_pos(ns_dn, t) != get_global_pos(get_nspace(coord), t) ) {
						printf("ERROR at mu = %i in direction %i\n", mu, dir);
						return;
					}
					coord.x = x;

					//mu = 2
					mu = 2;
					dir = 1;
					ns_up = get_neighbor(ns, mu);
					coord.y = (y + 1) % NSPACE;
					if( get_global_pos(ns_up, t) != get_global_pos(get_nspace(coord), t) ) {
						printf("ERROR at mu = %i in direction %i\n", mu, dir);
						return;
					}
					dir = -1;
					ns_dn = get_lower_neighbor(ns, mu);
					coord.y = (y - 1 + NSPACE) % NSPACE;
					if( get_global_pos(ns_dn, t) != get_global_pos(get_nspace(coord), t) ) {
						printf("ERROR at mu = %i in direction %i\n", mu, dir);
						return;
					}
					coord.y = y;

					//mu = 3
					mu = 3;
					dir = 1;
					ns_up = get_neighbor(ns, mu);
					coord.z = (z + 1) % NSPACE;
					if( get_global_pos(ns_up, t) != get_global_pos(get_nspace(coord), t) ) {
						printf("ERROR at mu = %i in direction %i\n", mu, dir);
						return;
					}
					dir = -1;
					ns_dn = get_lower_neighbor(ns, mu);
					coord.z = (z - 1 + NSPACE) % NSPACE;
					if( get_global_pos(ns_dn, t) != get_global_pos(get_nspace(coord), t) ) {
						printf("ERROR at mu = %i in direction %i\n", mu, dir);
						return;
					}
					coord.z = z;
				}
			}
		}
	}
}
