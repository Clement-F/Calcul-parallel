#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <filesystem>

using Mesh2DPart = std::vector<Mesh2D>;
// ========================== meme noeuds ===========================



std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega)
{  
    
    std::cout<<"\n instantiation of meshes \n";
    Mesh2DPart Meshes(4);   int sz = int(Omega.size());

    std::cout<<"\n creation of Q \n";
    CooMatrix<double> Q(sz,4);
    for(int j=0; j<sz; ++j){for(int k=0;k<4;k++){
        Q.push_back(j,k,0.);}}

    double j,k;
    auto nodes_= Omega.nodes();

    for(int i=0; i<4;i++){Meshes[i] = Mesh2D(nodes_);}

    long int coord=0;  
    double correction = 1.5;
    std::cout<<"\n begin second loop \n";
    for(int p =0;p<4;p++)    {
        j=p/2; k=p%2;
      for (auto el = Omega.begin(); el != Omega.end(); ++el)
      {
          auto element = *el;
          R3 Center = Ctr(element);
          if((Center[0]>(j/2 * correction) and Center[0]<((j+1)/2 * correction)) and (Center[1]>(k/2 * correction) and Center[1]<=((k+1)/2* correction))){Meshes[p].push_back(element); Q.push_back(coord,p,1); }
          coord++;
    }
    std::cout<<"\n mesh "<<p<<" has been sorted out \n";}
    
    std::cout<<"\n end looping \n";
    return std::make_tuple(Meshes, Q);
}



Mesh2D overlap(const Mesh2D& m_Omega, Mesh2D m_Gamma)
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
          not_in =false;    
        } 
      }     
      i++;
    }
    not_in=true; 
// -------------------------------------- 
  }
  // std::cout<<m_new_Gamma.size();
  return m_new_Gamma;
}


std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega, const std::size_t& nl)
{
  
    std::cout<<"\n instantiation of meshes \n";
    Mesh2DPart Meshes(4);   int sz = int(Omega.size());

    std::cout<<"\n creation of R \n";
    CooMatrix<double> R(sz,4);
    for(int j=0; j<sz; ++j){for(int k=0;k<4;k++){
        R.push_back(j,k,0.);}}

    double j,k;
    auto nodes_= Omega.nodes();

    for(int i=0; i<4;i++){Meshes[i] = Mesh2D(nodes_);}

    long int coord=0;  long int decalage=1;
    double correction = 1.5;
    std::cout<<"\n begin second loop \n";
    for(int p =0;p<4;p++)    {
        j=p/2; k=p%2;
      for (auto el = Omega.begin(); el != Omega.end(); ++el)
      {
          auto element = *el;
          R3 Center = Ctr(element);
          if((Center[0]>(j/2 * correction) and Center[0]<((j+1)/2 * correction)) and (Center[1]>(k/2 * correction) and Center[1]<=((k+1)/2* correction))){Meshes[p].push_back(element); R.push_back(coord,p,coord+decalage); }
          coord++;
    }decalage += coord; coord =0;
    std::cout<<"\n mesh "<<p<<" has been sorted out \n";}

    std::vector<Nodes> noeuds_partition(4);

    
    // tout les points sur lesquels se trouve les elements du maillage p
    for(int i =0;i<int(nodes_.size());i++)
    {
        R3 Point = nodes_[i];
        if(Point[0]<0.5)   {  if (Point[1]<0.5)     {noeuds_partition[0].push_back(Point);}
                                else                {noeuds_partition[1].push_back(Point);}}
        else                {  if (Point[1]<0.5)    {noeuds_partition[2].push_back(Point);}
                                else                {noeuds_partition[3].push_back(Point);}}
    }
    std::cout<<"\n ------------------------ \n";
    std::cout<<" le maillage possede "<<noeuds_partition[0].size()<<" noeuds et "<<Meshes[0].size()<<" elements ";
    std::cout<<"\n adding the "<<1<<" overlap of the "<<1<<"'s mesh  \n";
    Meshes[0]= overlap(Omega,Meshes[0]);

    std::cout<<" le maillage possede "<<noeuds_partition[0].size()<<" noeuds et "<<Meshes[0].size()<<" elements ";

    std::cout<<"\n adding the "<<2<<" overlap of the "<<1<<"'s mesh  \n";
    Meshes[0]= overlap(Omega,Meshes[0]);

    std::cout<<" le maillage possede "<<noeuds_partition[0].size()<<" noeuds et "<<Meshes[0].size()<<" elements ";

    std::cout<<"\n adding the "<<3<<" overlap of the "<<1<<"'s mesh  \n";
    Meshes[0]= overlap(Omega,Meshes[0]);

    std::cout<<" le maillage possede "<<noeuds_partition[0].size()<<" noeuds et "<<Meshes[0].size()<<" elements ";
    std::cout<<"\n ------------------------ \n";
  //   for(int k=0; k<int(nl);k++){

  //     for(int p =0; p<4;p++){ // à chaque maillage p 

  //     // ajoute tout les elements qui possède un point en commun avec noeuds_partition
  //     std::cout<<"\n adding the "<<k+1<<" overlap of the "<<p+1<<"'s mesh  \n";

  //     overlap(Omega,Meshes[p]);

  //     }
  // }



    std::cout<<"\n end looping \n";
    return std::make_tuple(Meshes, R);
}


std::pair< Mesh2DPart,CooMatrix<double> >
Partition16(const Mesh2D& Omega)
{  
    


    std::vector<bool> sorted(Omega.nodes().size(), false);
    std::vector<bool> mesh_sorted(Omega.size(), false);
    long int el_sorted=0; 

    double correction = 2./3;
    std::cout<<"\n instantiation of meshes \n";
    Mesh2DPart Meshes(16);   int sz = int(Omega.size());

    std::cout<<"\n creation of Q \n";
    CooMatrix<double> Q(sz,16);
    for(int j=0; j<sz; ++j){for(int k=0;k<16;k++){
        Q.push_back(j,k,0.);}}

    double j,k;
    std::cout<<"\n begin first loop \n";

    auto nodes_= Omega.nodes();
    for(int i=0; i<16;i++){Meshes[i] = Mesh2D(nodes_);}


    long int coord=0; 
    std::cout<<"\n begin second loop \n";
    for(int p =0;p<16;p++)    {
        j=p/4; k=p%4;
      for (auto el = Omega.begin(); el != Omega.end(); ++el)
      {
        // if(not mesh_sorted[coord]){
          auto element = *el; 
          R3 Center = Ctr(element);
        //   std::cout<<p<<'\t'<<Center<<'\t'<<j<<','<<k<<'\t';
          if(Center[0]*correction>j/4 and Center[0]*correction<(j+1)/4){if(Center[1]*correction>k/4 and Center[1]*correction<(k+1)/4){Meshes[p].push_back(element); Q.push_back(coord,p,1);}}                    
        // std::cout<<"\n-----------------\n";
        // }
        coord++;
    }
    el_sorted = std::count(mesh_sorted.begin(), mesh_sorted.end(), true);
    std::cout<<"\n mesh "<<p<<" has been sorted out \n";std::cout<<el_sorted<<" element have been sorted so far out of "<<sz<<'\n';
}
    
    std::cout<<"\n end looping \n";
    return std::make_tuple(Meshes, Q);
}

void
Plot(const std::vector<Mesh2D>& Sigma,const std::string& filename)
{
  std::filesystem::path filename_ = std::filesystem::path(filename);
  std::vector<std::string> tag = {"Vertices", "Edges", "Triangles", "Tetrahedra"};

  std::cout<<"\n create file\n";
  // Ouverture fichier
  filename_.replace_extension(".mesh");
  std::ofstream f;
  f.open(filename_.c_str());

  // Preambule
  f << "MeshVersionFormatted 3\n\n";
  f << "Dimension\n3\n\n";

    std::cout<<"\n view sizes \n";


  auto nodes_ = Sigma[0].nodes();
  const auto& v0 = nodes_[0];
  int el_size=0; 
  for(int i=0;i< int(Sigma.size());i++){el_size += Sigma[i].size();} 

  std::cout<<"\n begin first looping \n";

    
  // Section noeuds
  f << "Vertices\n";
  f << nodes_.size() << "\n";
  int k = 0;
//   auto v = Sigma[i].nodes();
  for(const auto& x:nodes_){f << x << "\t"<<1<<"\n"; k++;}
  
  f << "\n";

  
    std::cout<<"\n begin second looping \n";
    
  // Section elements
  f << tag[2]   << "\n";
  f << el_size << "\n";
  
  for(int i=0;i<int(Sigma.size());i++)
  {
  std::cout<<"\n the elements from mesh "<<i<<" have been plotted \n";
  for(const auto& e:Sigma[i]){
    for(std::size_t j=0; j<2+1; ++j){
      f << 1+int(&e[j]-&v0) << "\t";}
    f <<i+1 <<'\n';}
    }
  
  // Fermeture
  f << "\nEnd";
  f.close();
  
}


using FeSpace2D = FeSpace<2>;
using FeSpace2DxCoo = std::pair<FeSpace2D,CooMatrix<double>>;

FeSpace2DxCoo
Restrict(const FeSpace2D& Vh,const Mesh2D& Gamma,const std::vector<std::size_t>& tbl);

