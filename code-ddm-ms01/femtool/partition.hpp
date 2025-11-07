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
    Mesh2DPart Meshes(4);   int sz = int(Omega.size());

    std::cout<<"\n creation of Q \n";
    CooMatrix<double> Q(sz,4);
    for(int j=0; j<sz; ++j){for(int k=0;k<4;k++){
        Q.push_back(j,k,0.);}}

    int j =0;
    std::cout<<"\n begin first loop \n";
    auto nodes= Omega.nodes(); std::vector<Nodes> Nodes_part(4);
    for(int i =0;i<int(Omega.nodes().size());i++)
    {
        R3 Point = nodes[i];
        if(Point[0]<0.5)   {  if (Point[1]<0.5)     {Nodes_part[0].push_back(Point);}
                                else                {Nodes_part[1].push_back(Point);}}
        else                {  if (Point[1]<0.5)    {Nodes_part[2].push_back(Point);}
                                else                {Nodes_part[3].push_back(Point);}}
    }

    for(int i=0; i<4;i++){Meshes[i] = Mesh2D(Nodes_part[i]);}


    std::cout<<"\n begin second loop \n";
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


int 
find_index(const SmallVector<double,3>& ei, const Nodes& Nodes_) // brutefoooooooorce
{
    for(int i=0; i<int(Nodes_.size()); i++){
        if(Close(Nodes_[i],ei)){std::cout<<Nodes_[i]<<','<<ei<<'\n';   return i; };
    }
    return 0;
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

  int vert_size=0; int el_size=0;
  std::vector<Nodes> N(Sigma.size()); 
  for(int i=0;i< int(Sigma.size());i++){
    N[i] = (Sigma[i].nodes());  vert_size += N[i].size();  el_size += Sigma[i].size();} 

    std::cout<<"\n begin first looping \n";

    
  // Section noeuds
  f << "Vertices\n";
  f << vert_size << "\n";
  int k = 0;
  for(int i=0;i<int(Sigma.size());i++)
  {
//   auto v = Sigma[i].nodes();
  for(const auto& x:N[i]){

    // std::cout<<k<<" vectex, coord : "<<x<<'\n';
    f << x << "\t"<<i+1<<"\n"; k++;}
  }
  f << "\n";

  
    std::cout<<"\n begin second looping \n";
    
  // Section elements
  f << tag[2]   << "\n";
  f << el_size << "\n";
  
  int p, index, temp, decalage;

  for(int i=0;i<int(Sigma.size());i++)
  {

  p=0;  
  for(const auto& e:Sigma[i]){
    for(std::size_t j=0; j<2+1; ++j){

        index = 0;      temp=0; 
        decalage = 0;   k=0;
        while(index ==0 && k<int(Sigma.size()))
            { 
            std::cout<<p<<'\n';
            temp = find_index(e[j],N[k]); 
            if(temp!=0 or p==0){index = temp + decalage;std::cout<<temp<<'\n';}; 
            decalage += N[k].size(); 
            k++;
            }

        f << 1+ index << "\t";}
    f <<i+1<< "\n"; p++;}
    }
  // Fermeture
  f << "\nEnd";
  f.close();
  
}