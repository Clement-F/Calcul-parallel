#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <filesystem>

using Mesh2DPart = std::vector<Mesh2D>;



std::vector<std::size_t>
Restrict(const CooMatrix<double> R, const unsigned int p){
  size_t size = NbRow(R);
  auto new_R =MakeDense(R);
  std::vector<std::size_t> tbl; tbl.reserve(size);
  int i=0;
  for(int j=0;j<size;j++){if(new_R[j,p] !=0){tbl[i]=j; i++;}}
  return tbl;
}

//// ============================================================================
//// ========================== partition TP2 ===================================
//// ============================================================================


std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega)
{  
    
  std::vector<bool> mesh_sorted(Omega.size(), false);
  long int el_sorted=0; 

  Mesh2DPart Meshes(4);   int sz = int(Omega.size());

  CooMatrix<double> Q(sz,4);

  double j,k;
  auto nodes_= Omega.nodes();

  for(int i=0; i<4;i++){Meshes[i] = Mesh2D(nodes_);}

  long int coord=0;  
  double correction = 1.5;
  for(int p =0;p<4;p++)    {
      j=p/2; k=p%2;
    for (auto el = Omega.begin(); el != Omega.end(); ++el)
    {
      if(not mesh_sorted[coord]){
        auto element = *el;
        R3 Center = Ctr(element);
        if((Center[0]>(j/2 * correction) and Center[0]<((j+1)/2 * correction)) and (Center[1]>(k/2 * correction) and Center[1]<=((k+1)/2* correction)))
        {
          Meshes[p].push_back(element); 
          Q.push_back(coord,p,1); 
          mesh_sorted[coord]=true;
        }}
        coord++;
  }coord =0;
  el_sorted = std::count(mesh_sorted.begin(), mesh_sorted.end(), true);
  }
  
  return std::make_tuple(Meshes, Q);
}

std::pair< Mesh2DPart,CooMatrix<double> >
Partition16(const Mesh2D& Omega)
{  
  std::vector<bool> mesh_sorted(Omega.size(), false);
  long int el_sorted=0; 

  double correction = 2./3;
  Mesh2DPart Meshes(16);   int sz = int(Omega.size());

  CooMatrix<double> Q(sz,16);

  double j,k;

  auto nodes_= Omega.nodes();
  for(int i=0; i<16;i++){Meshes[i] = Mesh2D(nodes_);}


  long int coord=0; 
  for(int p =0;p<16;p++)    {
      j=p/4; k=p%4;
    for (auto el = Omega.begin(); el != Omega.end(); ++el)
    {
      if(not mesh_sorted[coord]){
        auto element = *el; 
        R3 Center = Ctr(element);
        if(Center[0]*correction>j/4 and Center[0]*correction<(j+1)/4){if(Center[1]*correction>k/4 and Center[1]*correction<(k+1)/4)
          {
            Meshes[p].push_back(element); 
            Q.push_back(coord,p,1);
            mesh_sorted[coord]=true;
          }}                    
      }
      coord++;
  }coord =0;
  el_sorted = std::count(mesh_sorted.begin(), mesh_sorted.end(), true);
  }
  
  return std::make_tuple(Meshes, Q);
}

Mesh2D overlap(const Mesh2D& m_Omega, Mesh2D m_Gamma, std::vector<std::size_t>& tbl)
{

  int sz = int(m_Omega.nodes().size());

  auto nodes_ = m_Omega.nodes();
  const auto& v0 = nodes_[0];

  std::vector<int> v(sz);
  int surcharge=0;
  for(const auto& el:m_Gamma){
      surcharge ++;
      if(surcharge>= sz){auto it = unique(v.begin(), v.end());    v.erase(it, v.end());}
      for(int j=0;j<3;j++){v.push_back(int(&el[j]-&v0));}
  }

  auto it = unique(v.begin(), v.end());
  v.erase(it, v.end());


  Mesh2D m_new_Gamma(m_Omega.nodes());
  int i=0; bool not_in = true;
  for(const auto& el:m_Omega){
    i=0;
// --------------------------------------
    while(i<int(v.size()) and not_in){  
      for(int j=0;j<3;j++){         
        if(int(&el[j]-&v0) ==v[i]){        
          m_new_Gamma.push_back(el); 
          tbl[j]=i;
          not_in =false;    
        } 
      }     
      i++;
    }
    not_in=true; 
// -------------------------------------- 
  }
  return m_new_Gamma;
}

std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega, const std::size_t& nl)
{
  
  std::vector<bool> mesh_sorted(Omega.size(), false);
  long int el_sorted=0; 

  Mesh2DPart Meshes(4);   int sz = int(Omega.size());

  CooMatrix<double> R(sz,4);

  double j,k;
  auto nodes_= Omega.nodes();

  for(int i=0; i<4;i++){Meshes[i] = Mesh2D(nodes_);}

  long int coord_global=0; 
  long int coord_loc   =0;  

  double correction = 1.5;
  for(int p =0;p<4;p++)    {
      j=p/2; k=p%2;
    for (auto el = Omega.begin(); el != Omega.end(); ++el)
    {
      coord_global++;
      if(not mesh_sorted[coord_global]){
        auto element = *el;
        R3 Center = Ctr(element);
        if((Center[0]>(j/2 * correction) and Center[0]<((j+1)/2 * correction)) and (Center[1]>(k/2 * correction) and Center[1]<=((k+1)/2* correction)))
        {
          Meshes[p].push_back(element); 
          R.push_back(coord_global,p,coord_loc); coord_loc++;
          mesh_sorted[coord_global]=true;
        }}
        
  }coord_global =0; coord_loc=0;
  // el_sorted = std::count(mesh_sorted.begin(), mesh_sorted.end(), true);
  }

  for(int k=0; k<int(nl);k++){

    for(int p =0; p<4;p++){ // à chaque maillage p 


    auto tbl_temp = Restrict(R,p);
    Meshes[p] = overlap(Omega,Meshes[p],tbl_temp);
    for(int i=0;i<sz;i++){if(tbl_temp[i]!=0) R.push_back(i,p,tbl_temp[i]);}

    }
  }
  return std::make_tuple(Meshes, R);
}

std::pair< Mesh2DPart,CooMatrix<double> >
Partition16(const Mesh2D& Omega, const std::size_t& nl)
{  
  
  std::vector<bool> mesh_sorted(Omega.size(), false);
  long int el_sorted=0; 

  double correction = 2./3;
  Mesh2DPart Meshes(16);   int sz = int(Omega.size());

  CooMatrix<double> R(sz,16);

  double j,k;

  auto nodes_= Omega.nodes();
  for(int i=0; i<16;i++){Meshes[i] = Mesh2D(nodes_);}


  long int coord=0; 
  for(int p =0;p<16;p++)    {
      j=p/4; k=p%4;
    for (auto el = Omega.begin(); el != Omega.end(); ++el)
    {
      if(not mesh_sorted[coord]){
        auto element = *el; 
        R3 Center = Ctr(element);
        if(Center[0]*correction>j/4 and Center[0]*correction<(j+1)/4){if(Center[1]*correction>k/4 and Center[1]*correction<(k+1)/4)
          {
            Meshes[p].push_back(element); 
            R.push_back(coord,p,coord);
            mesh_sorted[coord]=true;
          }}                    
      }
      coord++;
  }coord =0;
  el_sorted = std::count(mesh_sorted.begin(), mesh_sorted.end(), true);
  }
  
    
  for(int k=0; k<int(nl);k++){

    for(int p =0; p<4;p++){ // à chaque maillage p 


    auto tbl_temp = Restrict(R,p);
    Meshes[p] = overlap(Omega,Meshes[p],tbl_temp);
    for(int i=0;i<sz;i++){if(tbl_temp[i]!=0) R.push_back(i,p,tbl_temp[i]);}

    }
  }
  
  return std::make_tuple(Meshes, R);
}

void
Plot(const std::vector<Mesh2D>& Sigma,const std::string& filename)
{
  std::filesystem::path filename_ = std::filesystem::path(filename);
  std::vector<std::string> tag = {"Vertices", "Edges", "Triangles", "Tetrahedra"};

  // Ouverture fichier
  filename_.replace_extension(".mesh");
  std::ofstream f;
  f.open(filename_.c_str());

  // Preambule
  f << "MeshVersionFormatted 3\n\n";
  f << "Dimension\n3\n\n";



  auto nodes_ = Sigma[0].nodes();
  const auto& v0 = nodes_[0];
  int el_size=0; 
  for(int i=0;i< int(Sigma.size());i++){el_size += Sigma[i].size();} 


    
  // Section noeuds
  f << "Vertices\n";
  f << nodes_.size() << "\n";
  int k = 0;
//   auto v = Sigma[i].nodes();
  for(const auto& x:nodes_){f << x << "\t"<<1<<"\n"; k++;}
  
  f << "\n";

  
    
  // Section elements
  f << tag[2]   << "\n";
  f << el_size << "\n";
  
  for(int i=0;i<int(Sigma.size());i++)
  {
  for(const auto& e:Sigma[i]){
    for(std::size_t j=0; j<2+1; ++j){
      f << 1+int(&e[j]-&v0) << "\t";}
    f <<i+1 <<'\n';}
    }
  
  // Fermeture
  f << "\nEnd";
  f.close();
  
}



//// ============================================================================
//// ===================    partition TP 3     ==================================
//// ============================================================================



using FeSpace2D = FeSpace<2>;
using FeSpace2DxCoo = std::pair<FeSpace2D,CooMatrix<double>>;


FeSpace2DxCoo
Restrict(const FeSpace2D& Vh,const Mesh2D& Gamma,const std::vector<std::size_t>& tbl){
  auto Uh = FeSpace(Gamma);

  CooMatrix<double> P(dim(Uh),dim(Vh));
  for(int i=0;i<dim(Uh);i++){P.push_back(i,tbl[i],1);}

  return std::make_pair(Uh,P);
}


std::vector<FeSpace2DxCoo>
Partition4(const FeSpace2D& Vh, const std::size_t& nl)
{
  Mesh2D Omega = Vh.mesh();
  auto temp = Partition4(Omega,nl);

  auto R = temp.second;
  std::vector<FeSpace2DxCoo> Fe_x_R(4);
  for(int i=0;i<4;i++){

  auto Gamma = temp.first[i];
  auto R = temp.second;
  auto tbl = Restrict(R,0);

  Fe_x_R[i] = Restrict(Vh,Gamma,tbl);
  
  }

  return Fe_x_R;
}


std::vector<FeSpace2DxCoo>
Partition16(const FeSpace2D& Vh, const std::size_t& nl)
{
  Mesh2D Omega = Vh.mesh();
  auto temp = Partition16(Omega,nl);

  auto R = temp.second;
  std::vector<FeSpace2DxCoo> Fe_x_R(4);
  for(int i=0;i<16;i++){

  auto Gamma = temp.first[i];
  auto R = temp.second;
  auto tbl = Restrict(R,0);

  Fe_x_R[i] = Restrict(Vh,Gamma,tbl);
  
  }

  return Fe_x_R;
}