// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

//------------------------------------------//
//      FUNCTIONS INPUT PARSING             //
//------------------------------------------//
struct split_t
{
  enum empties_t { empties_ok, no_empties };
};

template <typename Container>
Container& split(Container& result,const typename Container::value_type& s,const typename Container::value_type& delimiters,split_t::empties_t empties = split_t::empties_ok)
{
  // splits the strings into individual fields
  // useful if you want to pass many parameters in the same string
  result.clear();
  size_t current;
  size_t next = -1;
  do
  {
    if (empties == split_t::no_empties)
    {
      next = s.find_first_not_of( delimiters, next + 1 );
      if (next == Container::value_type::npos) break;
      next -= 1;
    }
    current = next + 1;
    next = s.find_first_of( delimiters, current );
    result.push_back( s.substr( current, next - current ) );
  }
  while (next != Container::value_type::npos);
  return result;
}

void trim( std::string& s )
{
  // Remove leading and trailing whitespace
  static const char whitespace[] = " \n\t\v\r\f";
  s.erase( 0, s.find_first_not_of(whitespace) );
  s.erase( s.find_last_not_of(whitespace) + 1U );
}
//------------------------------------------//
//      end of FUNCTIONS INPUT PARSING      //
//------------------------------------------//


// function to check if file exists
inline bool fileExists(const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}



//invert matrix (for covariace calculation)
gsl_matrix *
invert_a_matrix(gsl_matrix *matrix,size_t size)
{
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}


float calculateFloodZ(UShort_t  *charge,
                      Crystal_t crystal)
{
  float FloodZ;
  float centralChargeOriginal;
  float centralSaturation;
  float centralPedestal;
  float division = 0.0;

  centralChargeOriginal = charge[crystal.detectorChannel];
  for(unsigned int iSat = 0; iSat < crystal.detectorSaturation.size(); iSat++)
  {
    if( crystal.detectorSaturation[iSat].digitizerChannel  == crystal.detectorChannel)
    {
      centralSaturation = crystal.detectorSaturation[iSat].saturation;
      centralPedestal = crystal.detectorSaturation[iSat].pedestal;
    }
  }
  float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );

  for (unsigned int iW = 0; iW < crystal.relevantForW.size(); iW++)
  {
    // std::cout << crystal.relevantForW[iW] << std::endl;
    float originalCh = charge[crystal.relevantForW[iW]];

    float saturationCh;
    float pedestalCorr;
    for(unsigned int iSat = 0; iSat < crystal.detectorSaturation.size(); iSat++)
    {
      if( crystal.detectorSaturation[iSat].digitizerChannel  == crystal.relevantForW[iW])
      {
        saturationCh = crystal.detectorSaturation[iSat].saturation;
        pedestalCorr = crystal.detectorSaturation[iSat].pedestal;
      }
    }
    // std::cout << originalCh << " "
    //           << saturationCh << " "
    //           << pedestalCorr << " "
    //           << std::endl;
    division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
  }

  FloodZ = centralChargeCorr / division;
  return FloodZ;
}
//calculate FloodZ...
