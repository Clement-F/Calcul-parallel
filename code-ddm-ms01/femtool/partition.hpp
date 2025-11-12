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

    int coord =0;
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
        std::cout<<'\t'<<Center<<'\t';
        if(Center[0]<0.75)   {  if (Center[1]<0.75)   {Meshes[0].push_back(element); Q.push_back(coord,0,1.);     }
                                else                {Meshes[1].push_back(element); Q.push_back(coord,1,1.);     }}
        else                {  if (Center[1]<0.75)   {Meshes[2].push_back(element); Q.push_back(coord,2,1.);     }
                                else                {Meshes[3].push_back(element); Q.push_back(coord,3,1.);     }}
        coord++;std::cout<<"\n-----------------\n";
    }
    
    std::cout<<"\n end looping \n";
    return std::make_tuple(Meshes, Q);
}

std::pair< Mesh2DPart,CooMatrix<double> >
Partition4_bis(const Mesh2D& Omega)
{  
    
    std::cout<<"\n instantiation of meshes \n";
    Mesh2DPart Meshes(4);   int sz = int(Omega.size());

    std::cout<<"\n creation of Q \n";
    CooMatrix<double> Q(sz,4);
    for(int j=0; j<sz; ++j){for(int k=0;k<4;k++){
        Q.push_back(j,k,0.);}}

    double j,k;
    std::cout<<"\n begin first loop \n";
    auto nodes= Omega.nodes(); std::vector<Nodes> Nodes_part(4);
    for(int p =0;p<4;p++)    {
        j=p/2; k=p%2;
        for(int i=0;i<int(Omega.nodes().size());i++)
        {
        R3 Point = nodes[i];
        if((Point[0]>j/2 and Point[0]<(j+1)/2)){if(Point[1]>k/2 and Point[1]<(k+1)/2){Nodes_part[p].push_back(Point);};}
        }
    }

    for(int i=0; i<4;i++){Meshes[i] = Mesh2D(Nodes_part[i]);}

    long int coord=0; 
    std::cout<<"\n begin second loop \n";
    for(int p =0;p<4;p++)    {
        j=p/2; k=p%2;
      for (auto el = Omega.begin(); el != Omega.end(); ++el)
      {
          Element<2> element = *el; 
          R3 Center = Ctr(element);
        //   std::cout<<p<<'\t'<<Center<<'\t'<<j<<','<<k<<'\t';
          if((Center[0]>j/2 and Center[0]<(j+1)/2) and (Center[1]>k/2 and Center[1]<(k+1)/2)){std::cout<<"passed";Meshes[p].push_back(element); Q.push_back(coord,p,coord+1); }
          coord++;std::cout<<"\n-----------------\n";
    }coord =0; 
    std::cout<<"\n mesh "<<p<<" has been sorted out \n";}
    
    std::cout<<"\n end looping \n";
    return std::make_tuple(Meshes, Q);
}

std::pair< Mesh2DPart,CooMatrix<double> >
Partition16(const Mesh2D& Omega)
{  
    


    std::vector<bool> sorted(Omega.nodes().size(), false);
    std::vector<bool> mesh_sorted(Omega.size(), false);
    long int el_sorted=0; long int node_sorted=0;

    double correction = 2./3;
    std::cout<<"\n instantiation of meshes \n";
    Mesh2DPart Meshes(16);   int sz = int(Omega.size());

    std::cout<<"\n creation of Q \n";
    CooMatrix<double> Q(sz,16);
    for(int j=0; j<sz; ++j){for(int k=0;k<16;k++){
        Q.push_back(j,k,0.);}}

    double j,k;
    std::cout<<"\n begin first loop \n";


    auto nodes= Omega.nodes(); std::vector<Nodes> Nodes_part(16);
    for(int p =0;p<16;p++)    {
        j=p/4; k=p%4;
        for(int i=0;i<int(Omega.nodes().size());i++)
        {
        // if(not sorted[i]){
        R3 Point = nodes[i];

                if((Point[0]>j/4 and Point[0]<(j+1)/4))         {if(Point[1]>k/4 and Point[1]<(k+1)/4){Nodes_part[p].push_back(Point);};    sorted[i] =true;}
        // else{   if(Point[0]==0)                                 {if(Point[1]>k/4 and Point[1]<(k+1)/4){Nodes_part[p].push_back(Point);};    sorted[i] =true;}
        // else{   if(Point[0]==1)                                 {if(Point[1]>k/4 and Point[1]<(k+1)/4){Nodes_part[p].push_back(Point);};    sorted[i] =true;}
        // else{   if((Point[0]>j/4 and Point[0]<(j+1)/4))         {if(Point[1]==0){Nodes_part[p].push_back(Point);};                          sorted[i] =true;}
        // else{   if((Point[0]>j/4 and Point[0]<(j+1)/4))         {if(Point[1]==1){Nodes_part[p].push_back(Point);};                          sorted[i] =true;}
        // }}}}
        // }
    }
        node_sorted = std::count(sorted.begin(), sorted.end(), true);
    std::cout<<"\n mesh "<<p<<" has been sorted out \n";std::cout<<node_sorted<<" nodes have been sorted so far out of "<<Omega.nodes().size()<<'\n';
    }

    for(int i=0; i<16;i++){Meshes[i] = Mesh2D(Nodes_part[i]);}

    long int coord=0; 
    std::cout<<"\n begin second loop \n";
    for(int p =0;p<16;p++)    {
        j=p/4; k=p%4;
      for (auto el = Omega.begin(); el != Omega.end(); ++el)
      {
        // if(not mesh_sorted[coord]){
          Element<2> element = *el; 
          R3 Center = Ctr(element);
        //   std::cout<<p<<'\t'<<Center<<'\t'<<j<<','<<k<<'\t';
          if(Center[0]*correction>j/4 and Center[0]*correction<(j+1)/4){if(Center[1]*correction>k/4 and Center[1]*correction<(k+1)/4){Meshes[p].push_back(element); Q.push_back(coord,p,coord+1);mesh_sorted[coord]= true;}}                    
        // std::cout<<"\n-----------------\n";
        // }
        coord++;
    }coord =0; 
    el_sorted = std::count(mesh_sorted.begin(), mesh_sorted.end(), true);
    std::cout<<"\n mesh "<<p<<" has been sorted out \n";std::cout<<el_sorted<<" element have been sorted so far out of "<<sz<<'\n';
}
    
    std::cout<<"\n end looping \n";
    return std::make_tuple(Meshes, Q);
}

// ========================== meme noeuds ===========================



std::pair< Mesh2DPart,CooMatrix<double> >
Partition4_ter(const Mesh2D& Omega)
{  
    
    std::cout<<"\n instantiation of meshes \n";
    Mesh2DPart Meshes(4);   int sz = int(Omega.size());

    std::cout<<"\n creation of Q \n";
    CooMatrix<double> Q(sz,4);
    for(int j=0; j<sz; ++j){for(int k=0;k<4;k++){
        Q.push_back(j,k,0.);}}

    double j,k;
    std::cout<<"\n begin first loop \n";
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
          if((Center[0]>(j/2 * correction) and Center[0]<((j+1)/2 * correction)) and (Center[1]>(k/2 * correction) and Center[1]<=((k+1)/2* correction))){Meshes[p].push_back(element); Q.push_back(coord,p,coord+1); }
          coord++;
    }coord =0; 
    std::cout<<"\n mesh "<<p<<" has been sorted out \n";}
    
    std::cout<<"\n end looping \n";
    return std::make_tuple(Meshes, Q);
}











int 
find_index(const SmallVector<double,3>& ei, const Nodes& Nodes_) // brutefoooooooorce
{
    for(int i=0; i<int(Nodes_.size()); i++){
        if(Close(Nodes_[i],ei)){return i; };
    }
    return -1;
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

    f << x << "\t"<<i+1<<"\n"; k++;}
  }
  f << "\n";

  
    std::cout<<"\n begin second looping \n";
    
  // Section elements
  f << tag[2]   << "\n";
  f << el_size << "\n";
  
  int index, temp, decalage;

  for(int i=0;i<int(Sigma.size());i++)
  {
  std::cout<<"\n the elements from mesh "<<i<<" have been plotted \n";
  for(const auto& e:Sigma[i]){
    for(std::size_t j=0; j<2+1; ++j){
        index = 0;      temp=0; 
        decalage = 0;   k=0;
        while(index ==0 && k<int(Sigma.size()))
            { 
            temp = find_index(e[j],N[k]); 
            if(temp!=-1){index = temp + decalage;}; 
            decalage += N[k].size(); 
            k++;
            }

        f << 1+ index << "\t";}
    f <<i+1<< "\n";}
    }
  // Fermeture
  f << "\nEnd";
  f.close();
  
}

void
Plot_bis(const std::vector<Mesh2D>& Sigma,const std::string& filename)
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


  auto nodes_ = Sigma[1].nodes();
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



