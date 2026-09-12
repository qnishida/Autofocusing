// Frozen loader reference from 51f60a7; do not update with production optimizations.
// Uses the public STATION representation, retaining the original loading/filter order.
namespace io_reference {
static std::vector<std::string> component_names(hid_t file_id, const std::string &station) {
  auto has = [&](const std::string &component) {
    return H5LTpath_valid(file_id, ("/" + station + "/" + component + "/sgram").c_str(), 1) > 0;
  };
  if (STATION::horizontal_only) {
    if (has("LE") && has("LN")) return {"LE", "LN"};
    if (has("E") && has("N")) return {"E", "N"};
  } else if (has("U") && has("E") && has("N")) {
    return {"U", "E", "N"};
  }
  return {};
}

void init_station(const std::string h5_file,std::vector<STATION> &sta0,const double rad0,const double rad1, hid_t fapl){//rad0 > rad1
  hid_t file_id = H5Fopen(h5_file.c_str(),H5F_ACC_RDONLY,fapl);//H5P_DEFAULT); //Open the file if it exists.
  if (file_id < 0) throw std::runtime_error("Cannot open HDF5: " + h5_file);
  hsize_t sta_num;
  H5Gget_num_objs(file_id,&sta_num);

  for(int i=0;i< (int)(sta_num);++i){
    char ctmp[100];
    H5Gget_objname_by_idx(file_id,(hsize_t)i, ctmp,sizeof(ctmp));
    std::string stnm = ctmp;
    
    if(H5Gget_objtype_by_idx(file_id,(hsize_t)i)==0 &&
       (!STATION::horizontal_only || !component_names(file_id, stnm).empty())){
      float stlo, stla,stel;
      H5LTget_attribute_float(file_id,("/"+stnm).c_str(),"stlo",&stlo);
      H5LTget_attribute_float(file_id,("/"+stnm).c_str(),"stla",&stla);
      H5LTget_attribute_float(file_id,("/"+stnm).c_str(),"stel",&stel);
      
      STATION sbuf;
      sbuf.set_station("Hi-net",(std::string) stnm,stla,stlo,stel);
      sta0.push_back(sbuf);
    }
  }
  H5Fclose(file_id);

  //Selection of stations and definition of station locations in Cartesian coordinate.
  //if(STATION::flag_ary==0){
  double arrlo=0.,arrla=0.;
  int iarr=0;
  //Initial location of the center of the array
  //double center_lon = 133., center_lat = 34.5;
  double center_lon = 138.5, center_lat = 36.5;
  for(std::vector<STATION>::iterator it = sta0.begin();it!=sta0.end();){
    const GeoDistance geo = calc_geodesic(center_lat, center_lon,
                                          it->print_lat(), it->print_lon());
    
    //cout << it->print_lat()<<" "<< it->print_lon()<<" "<< geo.dist_km<<" "<<rad0<<std::endl;
    if(geo.dist_km > rad0) it=sta0.erase(it);//Station selection
    else{
	if(geo.dist_km<rad1){
	  arrlo +=it->print_lon();
	  arrla +=it->print_lat();
	  iarr++;}
	it++;
    }
  }
  if(iarr == 0){
    std::ostringstream oss;
    oss << "No stations remain inside rad1=" << rad1
        << " km while initializing " << h5_file;
    throw std::runtime_error(oss.str());
  }
  STATION::lon_ary = arrlo/iarr;
  STATION::lat_ary = arrla/iarr;

  //cout<<STATION::lon_ary<<" "<<STATION::lat_ary<<std::endl;exit(0);
  for(int i=0;i < (int)sta0.size();i++){//Location of stations in Cartesian coordinate
    const GeoDistance geo = calc_geodesic(STATION::lat_ary, STATION::lon_ary,
                                          sta0[i].print_lat(), sta0[i].print_lon());
    
    double dx = geo.dist_km*sin(geo.az_deg/180.*M_PI);
    double dy = geo.dist_km*cos(geo.az_deg/180.*M_PI);
    sta0[i].set_loc(dx,dy,sqrt(dx*dx+dy*dy));
  }
}//exit(0);

int station_index(std::vector<STATION> &sta0,std::string sta,std::string net){
  int i=-1;

  //std::cerr << sta<<" "<<net<<std::endl;
  for(i=0;i<(int)sta0.size();i++){
    //std::cerr << sta0[i].print_sta()<<std::endl;
    //std::cerr << sta0[i].print_sta() <<"#"<<sta<<"#"<< (sta0[i].print_sta()==sta) <<std::endl;
    //std::cerr << sta0[i].print_net() <<"#"<<net<<"#"<< (sta0[i].print_net()==net) <<std::endl;
    if(sta0[i].print_sta()==sta&&sta0[i].print_net()==net) return(i);
  }
  return(-1);
}

int read_h5(std::vector<STATION> &sta0,std::string h5_file, hid_t fapl){
  int number=0;
  for (auto &station : sta0) station.clear_sac();
  int *ibuf = new int [STATION::npts];
  std::string stnm;
  SAC_data sac_buf;
  sac_buf.sgram = new float [STATION::npts];

  hid_t file_id = H5Fopen(h5_file.c_str(),H5F_ACC_RDONLY,fapl);//H5P_DEFAULT); //Open the file if it exists.
  if(file_id<0){
    delete[] ibuf;
    delete[] sac_buf.sgram;
    fprintf(stderr,"Cannot open file %s\n",h5_file.c_str());
    return(-1);
  }
  hsize_t sta_num;
  H5Gget_num_objs(file_id,&sta_num);

  for(int i=0;i< (int)(sta_num);i++){
    char cbuf[100];

    H5Gget_objname_by_idx(file_id,(hsize_t)i, cbuf,sizeof(cbuf));
    stnm = cbuf;
    const auto cmps = component_names(file_id, stnm);
    if(H5Gget_objtype_by_idx(file_id,(hsize_t)i)==0 && !cmps.empty()){
      float stlo, stla,stel;
      H5LTget_attribute_float(file_id,("/"+stnm).c_str(),"stlo",&stlo);
      H5LTget_attribute_float(file_id,("/"+stnm).c_str(),"stla",&stla);
      H5LTget_attribute_float(file_id,("/"+stnm).c_str(),"stel",&stel);

      int istnm = station_index(sta0,stnm,"Hi-net");//sta,net

      if(istnm !=-1){
	//components
	for(size_t j = 0;j<cmps.size();j++){
          int year=0,jday=0,hour=0,min=0,sec=0,msec=0,sr=0;
          float a0=0,cmpaz=0;
          const std::string group = "/"+stnm+"/"+cmps[j];
          const std::string time_group = group+"/time";
          if (H5LTget_attribute_float(file_id, group.c_str(), "sensitivity", &a0) < 0 ||
              H5LTget_attribute_float(file_id, group.c_str(), "cmpaz", &cmpaz) < 0 ||
              H5LTget_attribute_int(file_id, group.c_str(), "npts", &sac_buf.npts) < 0 ||
              H5LTget_attribute_int(file_id, group.c_str(), "sr", &sr) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "year", &year) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "jday", &jday) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "hour", &hour) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "min", &min) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "sec", &sec) < 0 ||
              H5LTget_attribute_int(file_id, time_group.c_str(), "msec", &msec) < 0 ||
              !std::isfinite(a0) || a0 == 0 || !std::isfinite(cmpaz)) {
            fprintf(stderr, "Invalid waveform metadata: %s\n", group.c_str());
            continue;
          }
          if (STATION::horizontal_only) {
            hsize_t unit_dims[1];
            H5T_class_t unit_class;
            size_t unit_size=0;
            if (H5LTget_attribute_info(file_id, group.c_str(), "unit", unit_dims,
                                      &unit_class, &unit_size) < 0 ||
                unit_class != H5T_STRING || unit_size > 32) {
              fprintf(stderr, "Missing/unsupported velocity unit: %s\n", group.c_str());
              continue;
            }
            char unit[33] = {};
            if (H5LTget_attribute_string(file_id, group.c_str(), "unit", unit) < 0 ||
                std::string(unit) != "nm/s") {
              fprintf(stderr, "Expected nm/s in %s\n", group.c_str());
              continue;
            }
          }
	  date d(date(year,1,1)+date_duration(jday-1));
	  sac_buf.ts = ptime(d,time_duration(hour,min,sec))+milliseconds(msec);
	  sac_buf.te = sac_buf.ts+milliseconds(sac_buf.npts*STATION::dt_msec);
	  sac_buf.Dt = 0;
	  sac_buf.cmpaz = cmpaz;
	  
	  int rank = 0;
          hsize_t dims[1] = {0};
          const std::string dataset = "/"+stnm+"/"+cmps[j]+"/sgram";
          if (H5LTget_dataset_ndims(file_id, dataset.c_str(), &rank) < 0 || rank != 1 ||
              H5LTget_dataset_info(file_id, dataset.c_str(), dims, nullptr, nullptr) < 0 ||
              dims[0] != static_cast<hsize_t>(sac_buf.npts)) {
            fprintf(stderr, "Invalid waveform extent: %s\n", dataset.c_str());
            continue;
          }
          if(1000 != sr*STATION::dt_msec || sac_buf.npts <= 0 || sac_buf.npts> STATION::npts){
	    fprintf(stderr,"Sampling rate (%d) or npts (%d) is wrong,\n",sr,sac_buf.npts);
	    sac_buf.npts = 0;
	  }
	  else{
	    if (H5LTread_dataset_int(file_id, dataset.c_str(), ibuf) < 0) continue;
	    for(int k=0;k<sac_buf.npts;k++) sac_buf.sgram[k]=ibuf[k]*a0*1E-9;
	    hp_filt(sac_buf.sgram,sac_buf.sgram,sac_buf.npts, 3E-2/(sr*1.));
            sta0[istnm].set_SAC_data(cmps[j].substr(cmps[j].size()-1),sac_buf);
	  }
	}
        if (STATION::horizontal_only && sta0[istnm].print_sta_num() == 2) {
          const double az_delta = std::fmod(sta0[istnm].print_sacE().cmpaz -
              sta0[istnm].print_sacN().cmpaz + 720., 360.);
          if (std::abs(az_delta - 90.) > 0.01) {
            fprintf(stderr, "Non-orthogonal LE/LN or E/N axes: %s\n", stnm.c_str());
            sta0[istnm].clear_sac();
          }
        }
        if (sta0[istnm].print_sta_num() == (STATION::horizontal_only ? 2 : 3)) ++number;
      }
    }
  }

  H5Fclose(file_id);
  delete[] ibuf;
  delete[] sac_buf.sgram;
  return(number);
}    
  
//init: location of the center of an array
//    : ibuf[86400*sr]?

}
