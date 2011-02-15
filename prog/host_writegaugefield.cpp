#include "host_writegaugefield.h"

hmc_error make_binary_data_single(ildg_gaugefield * array, char * out, const int array_size, const int num_entries){
  int length = sizeof(float);
  char * buf_tmp = new char[num_entries];
  int i, j;
  
  //save the array in a char-array, suppose it is one big array of array_size entries
  for(i = 0; i<array_size; i++){
    for(j = 0; j<length; j++){
      buf_tmp[i*length + j] = ( ((char*) &(array[0])) [i*length + j]  );
    }    
  }
  
  //suppose the buffer out has exactly the right size as given by array_size
  if(!ENDIAN){
    printf("\tThe ENDIANNESS of the system is little, bytes must be reversed\n");
    for (i=0; i<num_entries; i+=length){
      out[i]   = buf_tmp[i+3];
      out[i+1] = buf_tmp[i+2];
      out[i+2] = buf_tmp[i+1];
      out[i+3] = buf_tmp[i];
    }
  }
  else {
    printf("\tThe ENDIANNESS of the system is big, bytes must not be reversed\n");
    for (i=0; i<num_entries; i++){
      out[i] = buf_tmp[i];
    }
  }
  
  delete [] buf_tmp;
  return HMC_SUCCESS;
}

hmc_error make_binary_data_double(ildg_gaugefield * array, char * out, const int array_size, const int num_entries){
  int length = sizeof(double);
  char * buf_tmp = new char[num_entries];
  int i, j;

  //save the array in a char-array, suppose the array is one big array of array_size entries
  for(i = 0; i<array_size; i++){
    for(j = 0; j<length; j++){
      buf_tmp[i*length + j] = ( ((char*) &(array[0])) [i*length + j]  );
    }    
  }
  
  //suppose the buffer out has exactly the right size as given by array_size
  if(!ENDIAN){
    printf("\tThe ENDIANNESS of the system is little, bytes must be reversed\n");
    for (i=0; i<num_entries; i+=length){
      out[i]   = buf_tmp[i+7];
      out[i+1] = buf_tmp[i+6];
      out[i+2] = buf_tmp[i+5];
      out[i+3] = buf_tmp[i+4];
      out[i+4] = buf_tmp[i+3];
      out[i+5] = buf_tmp[i+2];
      out[i+6] = buf_tmp[i+1];
      out[i+7] = buf_tmp[i];
    }
  }
  else {
    printf("\tThe ENDIANNESS of the system is big, bytes must not be reversed\n");
    for (i=0; i<num_entries; i++){
      out[i] = buf_tmp[i];
    }
  }
  
  delete[] buf_tmp;
  return HMC_SUCCESS;
}


hmc_error write_gaugefield (
		ildg_gaugefield * array, int array_size, 
		int lx, int ly, int lz, int lt, int prec, int trajectorynr, hmc_float plaquettevalue, hmc_float beta, hmc_float kappa, hmc_float mu, hmc_float c2_rec, hmc_float epsilonbar, hmc_float mubar, 
		const char * hmc_version, const char * filename){

  printf("writing gaugefield to lime-file...\n");
  
  time_t current_time;
  FILE *outputfile;
  outputfile = fopen(filename, "w");
  int MB_flag;
  int ME_flag;
  n_uint64_t length_xlf_info=0, length_ildg_format=0, length_scidac_checksum = 0, length_ildg_binary_data = 0;
  LimeRecordHeader * header_ildg_format, *header_scidac_checksum, * header_ildg_binary_data, * header_xlf_info;

  //set values
  const char * field_out = "su3gauge";
  time_t rawtime;
  //here is the problem with current_time: rawtime has to be converted into a meaningful format
  current_time = rawtime;
  time ( &rawtime );
  const char * date = ctime (&rawtime);
  
  //get binary data
  //here it is assumed that the argument prec and sizeof(hmc_float) are the same!!
  int num_entries = sizeof(hmc_float)*array_size;

  char * binary_data = new char[num_entries];

  if(prec == 64){
    make_binary_data_double(array, binary_data, array_size, num_entries);
  }
  else if (prec == 32){
    make_binary_data_single(array, binary_data, array_size, num_entries);
  }
  else return HMC_STDERR;
  
  length_ildg_binary_data = num_entries;

  
  //write xlf-info to string, should look like this
  /*
  char xlf_info [] = "plaquette = 6.225960e-01\n trajectory nr = 67336\n beta = 6.000000, kappa = 0.177000, mu = 0.500000, c2_rec = 0.000000\n time = 1278542490\n hmcversion = 5.1.5\n mubar = 0.000000\n epsilonbar = 0.000000\n date = Thu Jul  8 00:41:30 2010\n";
  */
  
  //note: it is not clear what "time" is supposed to be
  char dummystring[1000];
  char xlf_info[1000];
  sprintf(xlf_info, "%s", "plaquette = ");
  sprintf(dummystring, "%f", plaquettevalue);
  strcat(xlf_info, dummystring);
  sprintf(dummystring, "%s", "\n trajectory nr = ");
  strcat(xlf_info, dummystring);
  strcat(xlf_info, "\n trajectory nr = ");
  sprintf(dummystring, "%i", trajectorynr);
  strcat(xlf_info, dummystring);  
  strcat(xlf_info, "\n beta = ");
  sprintf(dummystring, "%f", beta);
  strcat(xlf_info, dummystring);
  strcat(xlf_info, ", kappa = ");
  sprintf(dummystring, "%f", kappa);
  strcat(xlf_info, dummystring);  
  strcat(xlf_info, ", mu = ");
  sprintf(dummystring, "%f", mu);
  strcat(xlf_info, dummystring);  
  strcat(xlf_info, ", c2_rec = ");
  sprintf(dummystring, "%f", c2_rec);
  strcat(xlf_info, dummystring);  
  strcat(xlf_info, "\n time = ");
  sprintf(dummystring, "%i", (int) current_time);
  strcat(xlf_info, dummystring);
  strcat(xlf_info, "\n hmcversion = ");
  sprintf(dummystring, "%s", hmc_version);
  strcat(xlf_info, dummystring);  
  strcat(xlf_info, "\n mubar = ");
  sprintf(dummystring, "%f", mubar);
  strcat(xlf_info, dummystring);
  strcat(xlf_info, "\n epsilonbar = ");
  sprintf(dummystring, "%f", epsilonbar);
  strcat(xlf_info, dummystring);
  strcat(xlf_info, "\n date = ");
  sprintf(dummystring, "%s", date);
  strcat(xlf_info, dummystring);
 
  length_xlf_info = strlen(xlf_info);

  //write scidac checksum, this is stubb
  const char scidac_checksum [] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<scidacChecksum>\n  <version>1.0</version>\n  <suma>46b62a47</suma>\n  <sumb>1a24b4ac</sumb>\n</scidacChecksum>";
  length_scidac_checksum = strlen(scidac_checksum);
  
  //write ildg_format to string, should look like this:
  /*
  char ildg_format [] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n  <version>1.0</version>\n  <field>su3gauge</field>\n  <precision>64</precision>\n  <lx>4</lx>\n  <ly>4</ly>\n  <lz>4</lz>\n  <lt>4</lt>\n</ildgFormat>";
  */
  char dummystring2[1000];
  char ildg_format[1000];
  sprintf(ildg_format, "%s", "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n  <version>1.0</version>\n  <field>");
  strcat(ildg_format, field_out);
  strcat(ildg_format, "</field>\n  <precision>");
  sprintf(dummystring2, "%i", prec);
  strcat(ildg_format, dummystring2);
  strcat(ildg_format, "</precision>\n  <lx>");
  sprintf(dummystring2, "%i", lx);
  strcat(ildg_format, dummystring2);
  strcat(ildg_format, "</lx>\n  <ly>");
  sprintf(dummystring2, "%i", ly);
  strcat(ildg_format, dummystring2);
  strcat(ildg_format, "</ly>\n  <lz>");  
  sprintf(dummystring2, "%i", lz);
  strcat(ildg_format, dummystring2);
  strcat(ildg_format, "</lz>\n  <lt>");
  sprintf(dummystring2, "%i", lt);
  strcat(ildg_format, dummystring2);
  strcat(ildg_format, "</lt>\n</ildgFormat>");
  
  length_ildg_format = strlen(ildg_format);
  
  
  //write the lime file  
  MB_flag = 1; 
  
  LimeWriter *writer;
  writer = limeCreateWriter(outputfile);
  
  const char * types[] = {"xlf-info", "ildg-format", "ildg-binary-data", "scidac-checksum"};
  
  //xlf-info
  ME_flag = 1;
  header_xlf_info = limeCreateHeader(MB_flag, ME_flag, (char*) types[0], length_xlf_info);
  limeWriteRecordHeader(header_xlf_info, writer);
  limeDestroyHeader(header_xlf_info);
  limeWriteRecordData( xlf_info, &length_xlf_info, writer);
  printf("\txlf-info written\n");
 
  //ildg-format
  ME_flag = 2;
  header_ildg_format = limeCreateHeader(MB_flag, ME_flag, (char*) types[1], length_ildg_format);
  limeWriteRecordHeader(header_ildg_format, writer);
  limeDestroyHeader(header_ildg_format);
  limeWriteRecordData( ildg_format, &length_ildg_format, writer);  
  printf("\tildg-format written\n");

  //binary data
  ME_flag = 3;
  header_ildg_binary_data = limeCreateHeader(MB_flag, ME_flag, (char*) types[2], length_ildg_binary_data);
  limeWriteRecordHeader(header_ildg_binary_data, writer);
  limeDestroyHeader(header_ildg_binary_data);
  limeWriteRecordData(binary_data, &length_ildg_binary_data, writer);
  printf("\tildg_binary_data written\n");
  
  //scidac-checksum
  ME_flag = 4;
  header_scidac_checksum = limeCreateHeader(MB_flag, ME_flag, (char*) types[3], length_scidac_checksum);
  limeWriteRecordHeader(header_scidac_checksum, writer);
  limeDestroyHeader(header_scidac_checksum);
  limeWriteRecordData( (void*) scidac_checksum, &length_scidac_checksum, writer);
  printf("\tscidac-checksum written\n");
  
  //closing
  fclose(outputfile);
  limeDestroyWriter(writer);
  printf("\t%f MBytes were written to the lime file %s\n", (float) ( (float) (length_xlf_info+length_ildg_format + length_ildg_binary_data + length_scidac_checksum) /1024/1024 ),  filename);
  
  delete[] binary_data;
  return HMC_SUCCESS;
}
