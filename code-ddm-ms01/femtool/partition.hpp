#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <filesystem>

using Mesh2DPart = std::vector<Mesh2D>;

std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega)
{  
    
    std::cout<<"\n instantiation of meshes \n";
    Mesh2DPart Meshes(4);   std::size_t sz = Omega.size();

    std::cout<<"\n creation of Q \n";
    CooMatrix<double> Q(sz,4);
    for(std::size_t j=0; j<sz; ++j){for(int k=0;k<4;k++){
        Q.push_back(j,k,0.);}}

    int j =0;
    std::cout<<"\n begin looping \n";
    for (auto el = Omega.begin(); el != Omega.end(); ++el)
    {
        auto element = *el; 
        R3 Center = Ctr(element);
        // std::cout<<j<<" element, coord : "<<element<<'\n';
        if(Center[0]<0.5)   {  if (Center[1]<0.5)   {Meshes[0].push_back(element); Q.push_back(j,0,1.);     }
                                else                {Meshes[1].push_back(element); Q.push_back(j,1,1.);     }}
        else                {  if (Center[1]<0.5)   {Meshes[2].push_back(element); Q.push_back(j,2,1.);     }
                                else                {Meshes[3].push_back(element); Q.push_back(j,3,1.);     }}
        j++;
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

  int vert_size, el_size;
  std::vector<Nodes> N(Sigma.size()); 
  for(int i=0;i<Sigma.size();i++){N[i] = (Sigma[i].nodes());vert_size += N.size(); el_size += Sigma[i].size();} 

  const auto v0 = (N[0][0]);
  std::cout<<v0<<'\n';
    std::cout<<"\n begin first looping \n";

    
  // Section noeuds
  f << "Vertices\n";
  f << vert_size << "\n";
  int k = 0;
  for(int i=0;i<Sigma.size();i++)
  {
//   auto v = Sigma[i].nodes();
  for(const auto& x:N[i]){

    std::cout<<k<<" vectex, coord : "<<x<<'\n';
    f << x << "\t"<<i+1<<"\n"; k++;}
  f << "\n";
  }

  
    std::cout<<"\n begin second looping \n";
    
  // Section elements
  f << tag[2]   << "\n";
  f << el_size << "\n";
  for(int i=0;i<Sigma.size();i++)
  {
  for (auto el = Sigma[i].begin(); el != Sigma[i].end(); ++el){
    auto element = *el; 
    // std::cout<<k<<" element, coord : "<<element<<", "<<1+int(&element[0]-&v0)<< '\n';
    for(std::size_t j=0; j<2+1; ++j){
      f << 1+ int((&element[j]-&v0)) << "\t";}
    f <<i+1<< "\n";k++;}
    }
  // Fermeture
  f << "\nEnd";
  f.close();
  
}