
#ifdef _FERMIONS_

hmc_error opencl::init_fermion_variables(inputparameters* parameters, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	
	(*timer).reset();
		
	cout << "init solver variables..." << endl;
	int clerr = CL_SUCCESS; 
	int num_groups;
	if(local_work_size <= global_work_size) num_groups = global_work_size/local_work_size;
	else num_groups = 1;
	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	int eoprec_spinorfield_size = sizeof(hmc_complex)*EOPREC_SPINORFIELDSIZE;
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);
	int global_buf_size = complex_size*num_groups; 
	int global_buf_size_float = float_size*num_groups;
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
	hmc_complex kappa_complex = {(*parameters).get_kappa(), 0.};
	hmc_float tmp;
	
	cout << "\tinit spinorfields..." << endl;
  clmem_inout = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_source = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_source failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rn = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rn failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rhat = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rhat failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_v = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_p = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_p failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_s = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_s failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_t = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_t failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_aux = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //!!CP: eoprec fields, here one should insert a switch
	clmem_inout_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_source_even = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_source_odd = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clmem_rn_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rhat_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_v_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_p_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_s_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_t_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clmem_aux_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp_eoprec_1 = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clmem_tmp_eoprec_2 = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp_eoprec_3 = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  
  cout << "\tinit complex numbers..." << endl;
	clmem_rho = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rho failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rho_next = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rho_next failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_alpha = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_alpha failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_omega = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_omega failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_beta = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_beta failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp1 = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_tmp1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp2 = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_tmp2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_one = clCreateBuffer(context,CL_MEM_READ_ONLY,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_one failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_minusone = clCreateBuffer(context,CL_MEM_READ_ONLY,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_minusone failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_kappa_cmplx = clCreateBuffer(context,CL_MEM_READ_ONLY,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_kappa_cplx failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clmem_scalar_product_buf_glob = clCreateBuffer(context,CL_MEM_READ_WRITE,global_buf_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_scalar_product_buf_glob failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  

  cout << "\tinit float numbers..." << endl;
	clmem_kappa = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_kappa failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
	clmem_theta_fermion = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_theta_fermion failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_mu = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_mu failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_chem_pot_re = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_chem_pot_re failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_chem_pot_im = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_chem_pot_im failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_resid = clCreateBuffer(context,CL_MEM_READ_WRITE,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_resid failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
	clmem_trueresid = clCreateBuffer(context,CL_MEM_READ_WRITE,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_trueresid failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_global_squarenorm_buf_glob = clCreateBuffer(context,CL_MEM_READ_WRITE,global_buf_size_float,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_global_squarenorm_buf_glob failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  
  cout << "\twrite values to buffers..." << endl;
	tmp = (*parameters).get_kappa();
  clerr = clEnqueueWriteBuffer(queue,clmem_kappa,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_kappa failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  tmp = (*parameters).get_theta_fermion();
  clerr = clEnqueueWriteBuffer(queue,clmem_theta_fermion,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_theta_fermion failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  tmp = (*parameters).get_mu();
  clerr = clEnqueueWriteBuffer(queue,clmem_mu,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_mu failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  tmp = (*parameters).get_chem_pot_re();
  clerr = clEnqueueWriteBuffer(queue,clmem_chem_pot_re,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_chem_pot_re failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  tmp = (*parameters).get_chem_pot_im();
  clerr = clEnqueueWriteBuffer(queue,clmem_chem_pot_im,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_chem_pot_im failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

	clerr = clEnqueueWriteBuffer(queue,clmem_one,CL_TRUE,0,complex_size,&one,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_one failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  } 
  clerr = clEnqueueWriteBuffer(queue,clmem_minusone,CL_TRUE,0,complex_size,&minusone,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_minusone failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_kappa_cmplx,CL_TRUE,0,complex_size,&kappa_complex,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_kappa_cplx failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  cout << "\tinit fermion kernels..." << endl;
	M_diag = clCreateKernel(clprogram,"M_diag",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating M_diag kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	dslash = clCreateKernel(clprogram,"dslash",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating dslash kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	saxpy = clCreateKernel(clprogram,"saxpy",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating saxpy kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	saxsbypz = clCreateKernel(clprogram,"saxsbypz",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating saxsbypz kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	scalar_product = clCreateKernel(clprogram,"scalar_product",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating scalar_product kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	scalar_product_reduction = clCreateKernel(clprogram,"scalar_product_reduction",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating scalar_product_reduction kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	set_zero_spinorfield = clCreateKernel(clprogram,"set_zero_spinorfield",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating set_zero_spinorfield kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	global_squarenorm = clCreateKernel(clprogram,"global_squarenorm",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating global_squarenorm kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	global_squarenorm_reduction = clCreateKernel(clprogram,"global_squarenorm_reduction",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating global_squarenorm_reduction kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	ratio = clCreateKernel(clprogram,"ratio",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating ratio kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	product = clCreateKernel(clprogram,"product",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating product kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	convert_to_kappa_format = clCreateKernel(clprogram,"convert_to_kappa_format",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating convert_to_kappa_format kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	convert_from_kappa_format = clCreateKernel(clprogram,"convert_from_kappa_format",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating convert_from_kappa_format kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	convert_to_kappa_format_eoprec = clCreateKernel(clprogram,"convert_to_kappa_format_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating convert_to_kappa_format_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	convert_from_kappa_format_eoprec = clCreateKernel(clprogram,"convert_from_kappa_format_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating convert_from_kappa_format_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	create_point_source = clCreateKernel(clprogram,"create_point_source",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating create_point_source kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	M_sitediagonal = clCreateKernel(clprogram,"M_sitediagonal",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating M_sitediagonal kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	M_inverse_sitediagonal = clCreateKernel(clprogram,"M_inverse_sitediagonal",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating M_inverse_sitediagonal kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	dslash_eoprec = clCreateKernel(clprogram,"dslash_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating dslash_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	saxpy_eoprec = clCreateKernel(clprogram,"saxpy_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating saxpy_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	saxsbypz_eoprec = clCreateKernel(clprogram,"saxsbypz_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating saxsbypz_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	scalar_product_eoprec = clCreateKernel(clprogram,"scalar_product_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating scalar_product_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	set_zero_spinorfield_eoprec = clCreateKernel(clprogram,"set_zero_spinorfield_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating set_zero_spinorfield_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	global_squarenorm_eoprec = clCreateKernel(clprogram,"global_squarenorm_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating global_squarenorm_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	create_point_source_eoprec = clCreateKernel(clprogram,"create_point_source_eoprec",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating create_point_source_eoprec kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	
	
	//!!CP: this can be deleted later...
	cout << "\tset spinorfields to zero..." << endl;
	usetimer noop;
	set_zero_spinorfield_device(clmem_tmp, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_inout, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_source, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_rn, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_rhat, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_v, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_p, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_s, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_t, local_work_size, global_work_size, &noop);
	set_zero_spinorfield_device(clmem_aux, local_work_size, global_work_size, &noop);
 
  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::convert_to_kappa_format_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(convert_to_kappa_format,0,sizeof(cl_mem),&inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(convert_to_kappa_format,1,sizeof(cl_mem),&clmem_kappa); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clerr = clEnqueueNDRangeKernel(queue, convert_to_kappa_format,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue convert_to_kappa_format kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  
	(*timer).add();
	return HMC_SUCCESS;
}
	
hmc_error opencl::convert_from_kappa_format_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(convert_from_kappa_format,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(convert_from_kappa_format,1,sizeof(cl_mem),&out); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clerr = clSetKernelArg(convert_from_kappa_format,2,sizeof(cl_mem),&clmem_kappa); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clerr = clEnqueueNDRangeKernel(queue, convert_from_kappa_format,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue convert_from_kappa_format kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::convert_to_kappa_format_eoprec_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(convert_to_kappa_format_eoprec,0,sizeof(cl_mem),&inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(convert_to_kappa_format_eoprec,1,sizeof(cl_mem),&clmem_kappa); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clerr = clEnqueueNDRangeKernel(queue, convert_to_kappa_format_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue convert_to_kappa_format_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  
	(*timer).add();
	return HMC_SUCCESS;
}
	
hmc_error opencl::convert_from_kappa_format_eoprec_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(convert_from_kappa_format_eoprec,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(convert_from_kappa_format_eoprec,1,sizeof(cl_mem),&out); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clerr = clSetKernelArg(convert_from_kappa_format_eoprec,2,sizeof(cl_mem),&clmem_kappa); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clerr = clEnqueueNDRangeKernel(queue, convert_from_kappa_format_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue convert_from_kappa_format_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_spinorfield_to_device(hmc_spinor_field* host_spinorfield,  usetimer* timer){
// 	cout<<"Copy spinorfield to device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_inout,CL_TRUE,0,spinorfield_size,host_spinorfield,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_eoprec_spinorfield_to_device(hmc_eoprec_spinor_field* host_spinorfield,  usetimer* timer){
// 	cout<<"Copy spinorfield to device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*EOPREC_SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_inout_eoprec,CL_TRUE,0,spinorfield_size,host_spinorfield,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_source_to_device(hmc_spinor_field* host_source,  usetimer* timer){
// 	cout<<"Copy source to device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_source,CL_TRUE,0,spinorfield_size,host_source,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_eoprec_source_to_device(hmc_eoprec_spinor_field* host_source1, hmc_eoprec_spinor_field* host_source2, usetimer* timer){
// 	cout<<"Copy source to device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*EOPREC_SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_source_even,CL_TRUE,0,spinorfield_size,host_source1,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_source_odd,CL_TRUE,0,spinorfield_size,host_source2,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::get_spinorfield_from_device(hmc_spinor_field* host_spinorfield, usetimer* timer){
//   cout<<"Get spinorfield from device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
  int clerr = clEnqueueReadBuffer(queue,clmem_inout,CL_TRUE,0,spinorfield_size,host_spinorfield,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::get_eoprec_spinorfield_from_device(hmc_eoprec_spinor_field* host_spinorfield, usetimer* timer){
//   cout<<"Get spinorfield from device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*EOPREC_SPINORFIELDSIZE;
  int clerr = clEnqueueReadBuffer(queue,clmem_inout_eoprec,CL_TRUE,0,spinorfield_size,host_spinorfield,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::copy_spinor_device(cl_mem in, cl_mem out, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	
	clerr = clEnqueueCopyBuffer(queue,in, out, 0,0, spinorfield_size ,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_eoprec_spinor_device(cl_mem in, cl_mem out, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	int spinorfield_size = sizeof(hmc_complex)*EOPREC_SPINORFIELDSIZE;
	
	clerr = clEnqueueCopyBuffer(queue,in, out, 0,0, spinorfield_size ,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_float_from_device(cl_mem in, hmc_float * out, usetimer* timer){
	(*timer).reset();
	
	int clerr = CL_SUCCESS;
	hmc_float tmp;
	clerr = clEnqueueReadBuffer(queue,in,CL_TRUE,0,sizeof(hmc_float),&tmp,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }
	(*out) = tmp;
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_complex_from_device(cl_mem in, hmc_complex * out, usetimer* timer){
	(*timer).reset();
	
	int clerr = CL_SUCCESS;
	hmc_complex tmp;
	clerr = clEnqueueReadBuffer(queue,in,CL_TRUE,0,sizeof(hmc_complex),&tmp,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }
	(*out) = tmp;
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::M_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer, usetimer* dslashtimer, usetimer* Mdiagtimer){
	(*timer).reset();
	(*Mdiagtimer).reset();
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(M_diag,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,2,sizeof(cl_mem),&clmem_kappa);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,3,sizeof(cl_mem),&clmem_mu);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,M_diag,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue M_diag kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  (*Mdiagtimer).add();
	
	(*dslashtimer).reset();
	clerr = clSetKernelArg(dslash,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,1,sizeof(cl_mem),&clmem_tmp);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,2,sizeof(cl_mem),&clmem_gaugefield);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,3,sizeof(cl_mem),&clmem_theta_fermion);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,4,sizeof(cl_mem),&clmem_chem_pot_re);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,5,sizeof(cl_mem),&clmem_chem_pot_im);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,dslash,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue dslash kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	(*dslashtimer).add();
	
	//!! perhaps this can go into an extra calling of saxpy
	clerr = clSetKernelArg(saxpy,0,sizeof(cl_mem),&clmem_tmp); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,2,sizeof(cl_mem),&clmem_kappa_cmplx);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,3,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,saxpy,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue saxpy kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error opencl::Aee_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer, usetimer * singletimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * latimer){
	(*timer).reset();
	
	int even = EVEN;
	int odd = ODD;
	const size_t ls = local_work_size;
	const size_t gs = global_work_size;
	
	dslash_eoprec_device(in, clmem_tmp_eoprec_1, odd, ls, gs, dslashtimer);
	M_inverse_sitediagonal_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, ls, gs, Mdiagtimer);
	dslash_eoprec_device(clmem_tmp_eoprec_2, out, even, ls, gs, dslashtimer);
	M_sitediagonal_device(in, clmem_tmp_eoprec_1, ls, gs, Mdiagtimer);
	
	copy_eoprec_spinor_device(out, clmem_tmp_eoprec_3, singletimer);
	
	saxpy_eoprec_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_1, clmem_one, out, ls, gs, latimer);

// 	int clerr =CL_SUCCESS;
// 	clerr = clSetKernelArg(saxpy_eoprec,0,sizeof(cl_mem),&clmem_tmp_eoprec_3); 
//   if(clerr!=CL_SUCCESS) {
//     cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clerr = clSetKernelArg(saxpy_eoprec,1,sizeof(cl_mem),&clmem_tmp_eoprec_1);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clerr = clSetKernelArg(saxpy_eoprec,2,sizeof(cl_mem),&clmem_one);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clerr = clSetKernelArg(saxpy_eoprec,3,sizeof(cl_mem),&out);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clerr = clEnqueueNDRangeKernel(queue,saxpy_eoprec,1,0,&gs,&ls,0,0,NULL);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"enqueue saxpy_eoprec kernel failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::dslash_eoprec_device(cl_mem in, cl_mem out, int evenodd, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;	
	
	int eo = evenodd;
	clerr = clSetKernelArg(dslash_eoprec,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,2,sizeof(cl_mem),&clmem_gaugefield);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,3,sizeof(cl_mem),&clmem_theta_fermion);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,4,sizeof(cl_mem),&clmem_chem_pot_re);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,5,sizeof(cl_mem),&clmem_chem_pot_im);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,6,sizeof(cl_mem),&clmem_kappa);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 6 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,7,sizeof(int),&eo);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 7 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  clerr = clEnqueueNDRangeKernel(queue,dslash_eoprec,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue dslash_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);	
	
	return HMC_SUCCESS;
}

hmc_error opencl::M_inverse_sitediagonal_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	
	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(M_inverse_sitediagonal,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_inverse_sitediagonal,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_inverse_sitediagonal,2,sizeof(cl_mem),&clmem_kappa);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_inverse_sitediagonal,3,sizeof(cl_mem),&clmem_mu);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,M_inverse_sitediagonal,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue M_inverse_sitediagonal kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	
	return HMC_SUCCESS;
}

hmc_error opencl::M_sitediagonal_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	
	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(M_sitediagonal,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_sitediagonal,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_sitediagonal,2,sizeof(cl_mem),&clmem_kappa);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_sitediagonal,3,sizeof(cl_mem),&clmem_mu);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,M_sitediagonal,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue M_sitediagonal kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	
	return HMC_SUCCESS;
}

hmc_error opencl::saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(saxpy,0,sizeof(cl_mem),&x); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,1,sizeof(cl_mem),&y);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,2,sizeof(cl_mem),&alpha);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,3,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clEnqueueNDRangeKernel(queue,saxpy,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue saxpy kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}



hmc_error opencl::saxpy_eoprec_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(saxpy_eoprec,0,sizeof(cl_mem),&x); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy_eoprec,1,sizeof(cl_mem),&y);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy_eoprec,2,sizeof(cl_mem),&alpha);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy_eoprec,3,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clEnqueueNDRangeKernel(queue,saxpy_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue saxpy_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}
	
hmc_error opencl::saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(saxsbypz,0,sizeof(cl_mem),&x); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,1,sizeof(cl_mem),&y);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,2,sizeof(cl_mem),&z);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,3,sizeof(cl_mem),&alpha);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,4,sizeof(cl_mem),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,5,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,saxsbypz,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue saxsbypz kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::saxsbypz_eoprec_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(saxsbypz_eoprec,0,sizeof(cl_mem),&x); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz_eoprec,1,sizeof(cl_mem),&y);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz_eoprec,2,sizeof(cl_mem),&z);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz_eoprec,3,sizeof(cl_mem),&alpha);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz_eoprec,4,sizeof(cl_mem),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz_eoprec,5,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,saxsbypz_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue saxsbypz_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(scalar_product,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product,1,sizeof(cl_mem),&b);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(scalar_product,2,sizeof(cl_mem),&clmem_scalar_product_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product,3,sizeof(hmc_complex)*local_work_size,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  clerr = clEnqueueNDRangeKernel(queue,scalar_product,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  clerr = clSetKernelArg(scalar_product_reduction,0,sizeof(cl_mem),&clmem_scalar_product_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product_reduction,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,scalar_product_reduction,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product_reduction kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_complex_to_scalar_product_eoprec_device(cl_mem a, cl_mem b, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(scalar_product_eoprec,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product_eoprec,1,sizeof(cl_mem),&b);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product_eoprec,2,sizeof(cl_mem),&clmem_scalar_product_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product_eoprec,3,sizeof(hmc_complex)*local_work_size,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  clerr = clEnqueueNDRangeKernel(queue,scalar_product_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  clerr = clSetKernelArg(scalar_product_reduction,0,sizeof(cl_mem),&clmem_scalar_product_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product_reduction,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,scalar_product_reduction,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product_reduction kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  
  (*timer).add();
	return HMC_SUCCESS;
}


hmc_error opencl::set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(ratio,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(ratio,1,sizeof(cl_mem),&b);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(ratio,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //CP:this needs only one kernel!!
	size_t one = 1;
  clerr = clEnqueueNDRangeKernel(queue,ratio,1,0,&one,&one,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue ratio kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(product,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(product,1,sizeof(cl_mem),&b);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(product,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //CP:this needs only one kernel!!
	size_t one = 1;
  clerr = clEnqueueNDRangeKernel(queue,product,1,0,&one,&one,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue product kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_float_to_global_squarenorm_device(cl_mem a, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	
	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(global_squarenorm,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm,1,sizeof(cl_mem),&clmem_global_squarenorm_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(global_squarenorm,2,sizeof(hmc_float)*local_work_size,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  clerr = clEnqueueNDRangeKernel(queue,global_squarenorm,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue global_squarenorm kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  clerr = clSetKernelArg(global_squarenorm_reduction,0,sizeof(cl_mem),&clmem_global_squarenorm_buf_glob);      
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(global_squarenorm_reduction,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,global_squarenorm_reduction,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product_reduction kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_float_to_global_squarenorm_eoprec_device(cl_mem a, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	
	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(global_squarenorm_eoprec,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm_eoprec,1,sizeof(cl_mem),&clmem_global_squarenorm_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(global_squarenorm_eoprec,2,sizeof(hmc_float)*local_work_size,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  clerr = clEnqueueNDRangeKernel(queue,global_squarenorm_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue global_squarenorm kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  clerr = clSetKernelArg(global_squarenorm_reduction,0,sizeof(cl_mem),&clmem_global_squarenorm_buf_glob);      
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(global_squarenorm_reduction,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,global_squarenorm_reduction,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product_reduction kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_zero_spinorfield_device(cl_mem x, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(set_zero_spinorfield,0,sizeof(cl_mem),&x);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clEnqueueNDRangeKernel(queue,set_zero_spinorfield,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue set_zero_spinorfield kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_zero_spinorfield_eoprec_device(cl_mem x, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(set_zero_spinorfield_eoprec,0,sizeof(cl_mem),&x);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clEnqueueNDRangeKernel(queue,set_zero_spinorfield_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue set_zero_spinorfield_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_complex_device(cl_mem in, cl_mem out, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	int complex_size = sizeof(hmc_complex);
	
	clerr = clEnqueueCopyBuffer(queue,in, out, 0,0, complex_size ,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
	}	
	(*timer).add();
	return HMC_SUCCESS;
}
	
hmc_error opencl::bicgstab_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax){

	//!!CP: here one has to be careful if local_work_size is a null-pointer
	size_t globalsize = global_work_size;
	size_t localsize = local_work_size;
	//CP: these have to be on the host
	hmc_float resid;
	hmc_float trueresid;

	for(int iter=0; iter<cgmax; iter++){
		if(iter%iter_refresh==0) {
			set_zero_spinorfield_device(clmem_v, localsize, globalsize, latimer); 
			set_zero_spinorfield_device(clmem_p, localsize, globalsize, latimer);
			
			M_device(clmem_inout, clmem_rn, localsize, globalsize, Mtimer, dslashtimer, Mdiagtimer);
			saxpy_device(clmem_rn, clmem_source, clmem_one, clmem_rn, localsize, globalsize, latimer);
			copy_spinor_device(clmem_rn, clmem_rhat, singletimer);

			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			copy_complex_device(clmem_one, clmem_omega, singletimer);
			copy_complex_device(clmem_one, clmem_rho, singletimer);
			
			//CP: calc initial residuum for output, this is not needed for the algorithm!!
// 			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
// 			copy_float_from_device(clmem_resid, &resid, copytimer);
// 			cout << "initial residuum is: " << resid << endl;
		}

		set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1, singletimer);
		copy_complex_device(clmem_rho_next, clmem_rho, singletimer);
		set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2, singletimer);
		set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);

		set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1, singletimer);
		set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2, singletimer);
		saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p, local_work_size, global_work_size, latimer);

		M_device(clmem_p,clmem_v, local_work_size, global_work_size, Mtimer, dslashtimer, Mdiagtimer);

		set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha, singletimer);
		
		saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s, local_work_size, global_work_size, latimer);
		
		M_device(clmem_s, clmem_t, local_work_size, global_work_size, Mtimer, dslashtimer, Mdiagtimer);

		set_complex_to_scalar_product_device(clmem_t,clmem_s, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		//!!CP: this can also be global_squarenorm
		set_complex_to_scalar_product_device(clmem_t,clmem_t, clmem_tmp2, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega, singletimer);

		saxpy_device(clmem_t, clmem_s, clmem_omega, clmem_rn, local_work_size, global_work_size, latimer);

		saxsbypz_device(clmem_p, clmem_s, clmem_inout, clmem_alpha, clmem_omega, clmem_inout, local_work_size, global_work_size, latimer); 
    
		set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
		copy_float_from_device(clmem_resid, &resid, copytimer);

		if(resid<epssquare) {	
			M_device(clmem_inout,clmem_aux,local_work_size, global_work_size, Mtimer, dslashtimer, Mdiagtimer);
			saxpy_device(clmem_aux, clmem_source, clmem_one, clmem_aux, local_work_size, global_work_size, latimer); 
			set_float_to_global_squarenorm_device(clmem_aux, clmem_trueresid, local_work_size, global_work_size, scalarprodtimer);
			copy_float_from_device(clmem_trueresid, &trueresid, copytimer);
// 			cout << "residuum:\t" << resid << "\ttrueresiduum:\t" << trueresid << endl;
			if(trueresid<epssquare)
				return HMC_SUCCESS;
		}
		else{
// 			cout << "residuum:\t" << resid << endl;
		}

	}

	return HMC_SUCCESS;
}

hmc_error opencl::bicgstab_eoprec_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax){
	//!!CP: here one has to be careful if local_work_size is a null-pointer
	size_t globalsize = global_work_size;
	size_t localsize = local_work_size;
	//CP: these have to be on the host
	hmc_float resid;
	hmc_float trueresid;

	for(int iter=0; iter<cgmax; iter++){
		if(iter%iter_refresh==0) {
			set_zero_spinorfield_eoprec_device(clmem_v_eoprec, localsize, globalsize, latimer); 
			set_zero_spinorfield_eoprec_device(clmem_p_eoprec, localsize, globalsize, latimer);
			
			Aee_device(clmem_inout_eoprec, clmem_rn_eoprec, localsize, globalsize, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);
			saxpy_eoprec_device(clmem_rn_eoprec, clmem_source_even, clmem_one, clmem_rn_eoprec, localsize, globalsize, latimer);
			copy_eoprec_spinor_device(clmem_rn_eoprec, clmem_rhat_eoprec, singletimer);
			
			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			copy_complex_device(clmem_one, clmem_omega, singletimer);
			copy_complex_device(clmem_one, clmem_rho, singletimer);
			
			//!!CP: calc initial residuum, this is not needed for the algorithm!!
// 			set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
// 			copy_float_from_device(clmem_resid, &resid, copytimer);
// 			cout << "initial residuum is: " << resid << endl;
		}

		set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_rn_eoprec, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1, singletimer);
		copy_complex_device(clmem_rho_next, clmem_rho, singletimer);
		set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2, singletimer);
		set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);

		set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1, singletimer);
		set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2, singletimer);
		saxsbypz_eoprec_device(clmem_p_eoprec, clmem_v_eoprec, clmem_rn_eoprec, clmem_beta, clmem_tmp2, clmem_p_eoprec, local_work_size, global_work_size, latimer);

		Aee_device(clmem_p_eoprec,clmem_v_eoprec, local_work_size, global_work_size, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);
			
		set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_v_eoprec, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha, singletimer);
	
		saxpy_eoprec_device(clmem_v_eoprec, clmem_rn_eoprec, clmem_alpha, clmem_s_eoprec, local_work_size, global_work_size, latimer);
		
		Aee_device(clmem_s_eoprec, clmem_t_eoprec, local_work_size, global_work_size, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);
		
		set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec,clmem_s_eoprec, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		//!!CP: can this also be global_squarenorm??
		set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec,clmem_t_eoprec, clmem_tmp2, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega, singletimer);

		saxpy_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_omega, clmem_rn_eoprec, local_work_size, global_work_size, latimer);

		saxsbypz_eoprec_device(clmem_p_eoprec, clmem_s_eoprec, clmem_inout_eoprec, clmem_alpha, clmem_omega, clmem_inout_eoprec, local_work_size, global_work_size, latimer); 
    
		set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
		copy_float_from_device(clmem_resid, &resid, copytimer);
	
		if(resid<epssquare) {	
			Aee_device(clmem_inout_eoprec,clmem_aux_eoprec,local_work_size, global_work_size, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);
			saxpy_eoprec_device(clmem_aux_eoprec, clmem_source_even, clmem_one, clmem_aux_eoprec, local_work_size, global_work_size, latimer); 
			set_float_to_global_squarenorm_eoprec_device(clmem_aux_eoprec, clmem_trueresid, local_work_size, global_work_size, scalarprodtimer);
			copy_float_from_device(clmem_trueresid, &trueresid, copytimer);
// 			cout << "residuum:\t" << resid << "\ttrueresiduum:\t" << trueresid << endl;
			if(trueresid<epssquare)
				return HMC_SUCCESS;
		}
		else{
// 			cout << "residuum:\t" << resid << endl;
		}

	}

	return HMC_SUCCESS;
}

hmc_error opencl::cg_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax){
	//!!CP: here one has to be careful if local_work_size is a null-pointer
	size_t globalsize = global_work_size;
	size_t localsize = local_work_size;
	//CP: these have to be on the host
	hmc_float resid;
	int iter;
  for(iter = 0; iter < cgmax; iter ++){  
    if(iter%iter_refresh == 0){
			M_device(clmem_inout, clmem_rn, localsize, globalsize, Mtimer, dslashtimer, Mdiagtimer);
			saxpy_device(clmem_rn, clmem_source, clmem_one, clmem_rn, localsize, globalsize, latimer);
      copy_spinor_device(clmem_rn, clmem_p, copytimer);
    }
    //alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		set_complex_to_scalar_product_device(clmem_p, clmem_p, clmem_omega, local_work_size, global_work_size, scalarprodtimer);
		//A pn --> v
		M_device(clmem_p,clmem_v, local_work_size, global_work_size, Mtimer, dslashtimer, Mdiagtimer);
		set_complex_to_scalar_product_device(clmem_p, clmem_v, clmem_rho, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_omega, clmem_rho, clmem_alpha, singletimer);
		
		set_complex_to_product_device(clmem_alpha, clmem_minusone, clmem_tmp1, singletimer);

		saxpy_device(clmem_inout, clmem_p, clmem_tmp1, clmem_inout, localsize, globalsize, latimer);
		//rn+1 -> rhat
		saxpy_device(clmem_rn, clmem_v, clmem_alpha, clmem_rhat, localsize, globalsize, latimer);
		
		set_float_to_global_squarenorm_device(clmem_rhat, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
		copy_float_from_device(clmem_resid, &resid, copytimer);

		if(resid<epssquare) {
			return HMC_SUCCESS;
		}
		else{
			//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
			set_complex_to_scalar_product_device(clmem_rhat, clmem_rhat, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
			set_complex_to_ratio_device(clmem_rho_next, clmem_omega, clmem_beta, singletimer);
			
			//pn+1 = rn+1 + beta*pn 
			set_complex_to_product_device(clmem_beta, clmem_minusone, clmem_tmp2, singletimer);
			saxpy_device(clmem_p, clmem_rhat, clmem_tmp2, clmem_p, localsize, globalsize, latimer);
		}
	}
		return HMC_SUCCESS;
}
	


	
hmc_error opencl::simple_correlator_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * solvertimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t ls, const size_t gs, int cgmax){

	//!!CP: this has to be global or so
	int use_eo = 1;
	
	cout << "calc simple_propagator on the device..." << endl;
	hmc_spinor_field phi[SPINORFIELDSIZE];
	
	hmc_spinor_field in[SPINORFIELDSIZE];
	if(!use_eo){
		init_spinorfield_cold(in);
		copy_spinorfield_to_device(in, copytimer);
	}
	else{
		//!!CP: this should be fine since only half the field is used but of course it is not nice...
		init_spinorfield_cold_eoprec(in);
		copy_eoprec_spinorfield_to_device(in, copytimer);		
	}

  //pseudo scalar, flavour multiplet
  hmc_complex correlator_ps[NSPACE];
  for(int z=0; z<NSPACE; z++) {
    correlator_ps[z].re = 0;
    correlator_ps[z].im = 0;
  }
	
  for(int k=0; k<NC*NSPIN; k++) {
		if(!use_eo){
			create_point_source_device(k,0,0,ls, gs, latimer);
			solver_device(phi, copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, cgmax);
		}
		else{
			create_point_source_eoprec_device(k,0,0,ls, gs, latimer, dslashtimer, Mdiagtimer);
			solver_eoprec_device(phi, copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, cgmax);
		}
    for(int timepos = 0; timepos<NTIME; timepos++) {
			for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
				for(int alpha = 0; alpha<NSPIN; alpha++) {
					for(int c = 0; c<NC; c++) {
					// int j = spinor_element(alpha,c);
					int n = spinor_field_element(alpha, c, spacepos, timepos);
					int z = get_spacecoord(spacepos, 3);
					hmc_complex tmp = phi[n];
					hmc_complex ctmp = complexconj(&tmp);
					hmc_complex incr = complexmult(&ctmp,&tmp);
					correlator_ps[z].re += incr.re;
					correlator_ps[z].im += incr.im;
		}}}}
  }

  printf("pseudo scalar correlator:\n");
  for(int z=0; z<NSPACE; z++) {
    printf("%d\t(%e,%e)\n",z,correlator_ps[z].re,correlator_ps[z].im);
  }
  return HMC_SUCCESS;
}

hmc_error opencl::solver_device(hmc_spinor_field* out, usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax){
	(*solvertimer).reset();
	convert_to_kappa_format_device(clmem_inout, ls, gs, latimer);
	bicgstab_device(copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer ,ls, gs, cgmax);
	convert_from_kappa_format_device(clmem_inout, clmem_inout, ls, gs, latimer);
	get_spinorfield_from_device(out, copytimer);
	(*solvertimer).add();
	
	return HMC_SUCCESS;
}


hmc_error opencl::solver_eoprec_device(hmc_spinor_field* out, usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax){
	(*solvertimer).reset();
	hmc_eoprec_spinor_field phi_even[EOPREC_SPINORFIELDSIZE];
	hmc_eoprec_spinor_field phi_odd[EOPREC_SPINORFIELDSIZE];
	
	//CP: even solution
	convert_to_kappa_format_eoprec_device(clmem_inout_eoprec, ls, gs, latimer);
	bicgstab_eoprec_device(copytimer, singletimer, Mtimer, scalarprodtimer, latimer ,dslashtimer, Mdiagtimer, ls, gs, cgmax);	
	convert_from_kappa_format_eoprec_device(clmem_inout_eoprec, clmem_inout_eoprec, ls, gs, latimer);
	get_eoprec_spinorfield_from_device(phi_even, copytimer);

	//P: odd solution
	//!!CP: this transformation is not really necessary!!
	//!!CP: perhaps one can save some variables used here
	convert_to_kappa_format_eoprec_device(clmem_inout_eoprec, ls, gs, latimer);
	dslash_eoprec_device(clmem_inout_eoprec, clmem_tmp_eoprec_3, ODD, ls, gs, dslashtimer);
	M_inverse_sitediagonal_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_1, ls, gs, Mdiagtimer);
	M_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_2, ls, gs, Mdiagtimer);
	saxpy_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, clmem_one, clmem_inout_eoprec,  ls, gs, latimer);
	convert_from_kappa_format_eoprec_device(clmem_inout_eoprec, clmem_inout_eoprec,  ls, gs, latimer);
	get_eoprec_spinorfield_from_device(phi_odd, copytimer);

	//CP: whole solution
	convert_from_eoprec(phi_even,phi_odd,out);
	
	(*solvertimer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::create_point_source_device(int i, int spacepos, int timepos, const size_t ls, const size_t gs, usetimer * latimer){
	
	set_zero_spinorfield_device(clmem_source, ls, gs, latimer);
	
	(*latimer).reset();
	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(create_point_source,0,sizeof(cl_mem),&clmem_source);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(create_point_source,1,sizeof(int),&i);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(create_point_source,2,sizeof(int),&spacepos);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(create_point_source,3,sizeof(int),&timepos);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(create_point_source,4,sizeof(cl_mem),&clmem_kappa);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,create_point_source,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue dslash kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clFinish(queue);
	
	(*latimer).add();
	return HMC_SUCCESS;
}
	

hmc_error opencl::create_point_source_eoprec_device(int i, int spacepos, int timepos, const size_t ls, const size_t gs, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer){
	
	set_zero_spinorfield_eoprec_device(clmem_source_even, ls, gs, latimer);
	set_zero_spinorfield_eoprec_device(clmem_source_odd, ls, gs, latimer);
	set_zero_spinorfield_eoprec_device(clmem_tmp_eoprec_2, ls, gs, latimer);
	set_zero_spinorfield_eoprec_device(clmem_tmp_eoprec_1, ls, gs, latimer);

	//!!CP: I do not know why this does not work!! The compiler says he does not know get_global_pos, which does not make sense
	int glob_pos = spacepos + VOLSPACE * timepos;// get_global_pos(spacepos, timepos);
	int n = get_n_eoprec(spacepos, timepos);
	int clerr = CL_SUCCESS;
	int evenodd = glob_pos%2;

	(*latimer).reset();
	//CP: this is different than the host code, where this is done implicitly when converting the normal source to even/odd
	if(evenodd == 1){
		clerr = clSetKernelArg(create_point_source_eoprec,0,sizeof(cl_mem),&clmem_source_odd);
		if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
		}
  }
	else {
		clerr = clSetKernelArg(create_point_source_eoprec,0,sizeof(cl_mem),&clmem_tmp_eoprec_2);
		if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
		}		
	}
	clerr = clSetKernelArg(create_point_source_eoprec,1,sizeof(int),&i);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(create_point_source_eoprec,2,sizeof(int),&n);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(create_point_source_eoprec,3,sizeof(cl_mem),&clmem_kappa);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,create_point_source_eoprec,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue create_point_source_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }		
	clFinish(queue);
	(*latimer).add();

	M_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_1, ls, gs, Mdiagtimer);
	dslash_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_3, EVEN, ls, gs, dslashtimer);
	
	saxpy_eoprec_device(clmem_tmp_eoprec_2, clmem_tmp_eoprec_3, clmem_one, clmem_source_even, ls, gs, latimer);

	return HMC_SUCCESS;
}
	
hmc_error opencl::finalize_fermions(){
	
  if(clFlush(queue)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clFinish(queue)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseKernel(M_diag)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(dslash)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(saxpy)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(saxsbypz)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(scalar_product)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(scalar_product_reduction)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(set_zero_spinorfield)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(global_squarenorm)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(global_squarenorm_reduction)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(ratio)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(product)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(convert_to_kappa_format)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(convert_from_kappa_format)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(convert_to_kappa_format_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(convert_from_kappa_format_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(create_point_source)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(M_sitediagonal)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(M_inverse_sitediagonal)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(dslash_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(saxpy_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(scalar_product_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(set_zero_spinorfield_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(global_squarenorm_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(create_point_source_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);

																															
  if(clReleaseMemObject(clmem_kappa)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_theta_fermion)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_mu)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_chem_pot_re)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_chem_pot_im)!=CL_SUCCESS) exit(HMC_OCLERROR);
	
	if(clReleaseMemObject(clmem_inout)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_source)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rn)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rhat)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_v)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_p)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_s)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_t)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_aux)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp)!=CL_SUCCESS) exit(HMC_OCLERROR);
	
	if(clReleaseMemObject(clmem_inout_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_source_even)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_source_odd)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rn_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rhat_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_v_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_p_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_s_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_t_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_aux_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp_eoprec_1)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp_eoprec_2)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp_eoprec_3)!=CL_SUCCESS) exit(HMC_OCLERROR);
	
	if(clReleaseMemObject(clmem_rho)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rho_next)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_alpha)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_omega)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_beta)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp1)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp2)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_one)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_minusone)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_kappa_cmplx)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_scalar_product_buf_glob)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_global_squarenorm_buf_glob)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_resid)!=CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_trueresid)!=CL_SUCCESS) exit(HMC_OCLERROR);
  
	return HMC_SUCCESS;
}
	
hmc_error opencl::perform_benchmark(int steps, int cgmax, const size_t ls, const size_t gs, usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * solvertimer, usetimer * dslashtimer, usetimer * Mdiagtimer){
	hmc_float resid;
	hmc_spinor_field phi[SPINORFIELDSIZE];
	for(int i = 0; i<steps; i++){
		if(!use_eo){
			convert_to_kappa_format_device(clmem_inout, ls, gs, latimer);
			set_zero_spinorfield_device(clmem_v, ls, gs, latimer);
			saxpy_device(clmem_rn, clmem_source, clmem_one, clmem_rn, ls, gs, latimer);
			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next, ls, gs, scalarprodtimer);
			set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1, singletimer);
			saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p, ls, gs, latimer);
			M_device(clmem_inout, clmem_rn, ls, gs, Mtimer, dslashtimer, Mdiagtimer);
			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, ls, gs, scalarprodtimer);
			copy_float_from_device(clmem_resid, &resid, copytimer);
			copy_spinor_device(clmem_rn, clmem_rhat, singletimer);
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);
			create_point_source_device(0,0,0,ls, gs, latimer);
			solver_device(phi, copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, cgmax);
			convert_from_kappa_format_device(clmem_inout, clmem_inout, ls, gs, latimer);
			get_spinorfield_from_device(phi, copytimer);
		}
		else{
			set_zero_spinorfield_eoprec_device(clmem_v_eoprec, ls, gs, latimer); 
			Aee_device(clmem_inout_eoprec, clmem_rn_eoprec, ls, gs, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);
			copy_eoprec_spinor_device(clmem_rn_eoprec, clmem_rhat_eoprec, singletimer);
			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2, singletimer);
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);
			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_rn_eoprec, clmem_rho_next, ls, gs, scalarprodtimer);
			saxsbypz_eoprec_device(clmem_p_eoprec, clmem_v_eoprec, clmem_rn_eoprec, clmem_beta, clmem_tmp2, clmem_p_eoprec, ls, gs, latimer);
			saxpy_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_omega, clmem_rn_eoprec, ls, gs, latimer);
			set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid, ls, gs, scalarprodtimer);
			copy_float_from_device(clmem_resid, &resid, copytimer);
			create_point_source_eoprec_device(0,0,0,ls, gs, latimer, dslashtimer, Mdiagtimer);
			solver_eoprec_device(phi, copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, cgmax);
			convert_to_kappa_format_eoprec_device(clmem_inout_eoprec, ls, gs, latimer);
			convert_from_kappa_format_eoprec_device(clmem_inout_eoprec, clmem_inout_eoprec, ls, gs, latimer);
			get_eoprec_spinorfield_from_device(phi, copytimer);
		}
	}
	
	
	return HMC_SUCCESS;
}

#endif