#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <filesystem>
#include <set>
#include <utility>

using Mesh2DPart = std::vector<Mesh2D>;

int 
find_index(const SmallVector<double,3>& ei, const Nodes& Nodes_)
{
    for(int i=0; i<int(Nodes_.size()); i++){
        if(Close(Nodes_[i],ei)){return i; };
    }
    return -1;
}

std::pair<Mesh2DPart, CooMatrix<double>> 
 Partition4(const Mesh2D& Omega) {

    Mesh2DPart Sigma(4);
    int ne = int(Omega.size());
    CooMatrix<double> Q(ne, 4);

    // Boucle sur les éléments
    for (int j=0; j<ne; ++j) {
        auto elem = Omega[j];
        R3 centre = Ctr(elem);

        double x = centre[0];
        double y = centre[1];
        
        for (int p=0; p<4; ++p) {
            double j_idx = p/2;
            double k_idx = p % 2;

            double x_min =  (j_idx) / 2.0;
            double x_max =  (j_idx + 1)/ 2.0;
            double y_min = (k_idx) / 2.0;
            double y_max = (k_idx + 1) / 2.0;
            if (x > x_min && x < x_max && y > y_min && y < y_max) {
                Sigma[p].push_back(elem);
                Q.push_back(j, p, 1.0);
                break; // on a le sous-domaine, pas besoin de regarder les autres
            }
        }
    }

    // Reconstruire chaque sous-maillage avec une liste de noeuds locaux
    for (int p=0; p<4; ++p) {
        Nodes noeuds_submesh;
        for (const auto& elem : Sigma[p]) {
            for (int v=0; v<3; ++v) {
                const auto& noeud = elem[v];
                if (find_index(noeud, noeuds_submesh) == -1) {
                    noeuds_submesh.push_back(noeud);
                }
            }
        }
        Mesh2D reconstruction(noeuds_submesh);
        for (const auto& elem : Sigma[p]) {
            reconstruction.push_back(elem);
        }
        Sigma[p] = std::move(reconstruction);
    }
    return make_pair(Sigma, Q);
}

std::pair<Mesh2DPart, CooMatrix<double>>
 Partition4(const Mesh2D& Omega, const std::size_t& nl) {

    auto [Sigma, Q] = Partition4(Omega);
    int ne = Omega.size();

    // Construire les connections entre triangles (voisins, pour rajouter nl couches après)
    std::vector<std::vector<int>> voisins(ne);
    for (int i=0; i<ne; ++i) {
        for (int j=i+1; j<ne; ++j) {
            auto& elem_i = Omega[i];
            auto& elem_j = Omega[j];
            int connexions = 0;
            for (int a=0; a<3; ++a) {
                for (int b=0; b<3; ++b) {
                    if (Close(elem_i[a], elem_j[b])) {
                        connexions++;
                    }
                }
            }
            if (connexions >= 2) {
                voisins[i].push_back(j);
                voisins[j].push_back(i);
            }
        }
    }

    // Ajouter nl couches
    std::vector<std::set<int>> elements_Gamma(4);
    const auto& Qdata = GetData(Q);
    for (const auto& [j, pcol, val] : Qdata) {
        if (val != 0.0) {
            elements_Gamma[static_cast<int>(pcol)].insert(j);
        }
    }
    for (int p=0; p<4; ++p) {
        for (std::size_t couche=0; couche<nl; ++couche) {
            std::set<int> new_couche;
            for (auto j : elements_Gamma[p]) {
                for (auto v : voisins[j]) {
                    new_couche.insert(v);
                }
            }
        elements_Gamma[p].insert(new_couche.begin(), new_couche.end());
        }
    }

    // Construire Gamma et R
    Mesh2DPart Gamma(4);
    CooMatrix<double> R(ne, 4);
    for (int p=0; p<4; ++p) {
        int index_local = 1; // numérotation pour le sous domaine
        for (int j=0; j<ne; ++j) {
            if (elements_Gamma[p].count(j)) {
                Gamma[p].push_back(Omega[j]);
                R.push_back(j, p, index_local++);
            }
        }
    }

    // Reconstruire chaque sous-maillage avec une liste de noeuds locaux
    for (int p=0; p<4; ++p) {
        Nodes noeuds_submesh;
        for (const auto& elem : Gamma[p]) {
            for (int v=0; v<3; ++v) {
                const auto& arete = elem[v];
                if (find_index(arete, noeuds_submesh) == -1) {
                    noeuds_submesh.push_back(arete);
                }
            }
        }
        Mesh2D reconstruction(noeuds_submesh);
        for (const auto& elem : Gamma[p]) {
            reconstruction.push_back(elem);
        }
        Gamma[p] = std::move(reconstruction);
    }
    return make_pair(Gamma, R);
}

std::pair<Mesh2DPart, CooMatrix<double>> 
 Partition16(const Mesh2D& Omega) {

    int ne = int(Omega.size());
    Mesh2DPart Sigma(16);
    CooMatrix<double> Q(ne, 16);


    for (int j=0; j<ne; ++j) {
        auto elem = Omega[j];
        R3 centre = Ctr(elem);

        double x = centre[0];
        double y = centre[1];

        for (int p=0; p<16; ++p) {
            double j_idx = p/4;
            double k_idx = p % 4;

            double x_min = (j_idx) / 4;
            double x_max = (j_idx + 1) / 4;
            double y_min = (k_idx) / 4;
            double y_max = (k_idx + 1) / 4;
            if (x > x_min && x < x_max && y > y_min && y < y_max) {
                Sigma[p].push_back(elem);
                Q.push_back(j, p, 1.0);
                break;
            }
        } 
    }
    for (int p=0; p<16; ++p) {
        Nodes noeuds_submesh;
        for (const auto& elem : Sigma[p]) {
            for (int v=0; v<3; ++v) {
                const auto& noeud = elem[v];
                if (find_index(noeud, noeuds_submesh) == -1) {
                    noeuds_submesh.push_back(noeud);
                }
            }
        }
        Mesh2D reconstruction(noeuds_submesh);
        for (const auto& elem : Sigma[p]) {
            reconstruction.push_back(elem);
        }
        Sigma[p] = std::move(reconstruction);
    }
    return make_pair(Sigma, Q);
 }

std::pair<Mesh2DPart, CooMatrix<double>> 
 Partition16(const Mesh2D& Omega, const std::size_t& nl) {

    auto [Sigma, Q] = Partition16(Omega);
    int ne = Omega.size();

    // Construire les connexions entre triangles
    std::vector<std::vector<int>> voisins(ne);
    for (int i=0; i<ne; ++i) {
        for (int j=i+1; j<ne; ++j) {
            auto& elem_i = Omega[i];
            auto& elem_j = Omega[j];
            int connexions = 0;
            for (int a=0; a<3; ++a) {
                for (int b=0; b<3; ++b) {
                    if (Close(elem_i[a], elem_j[b])) {
                        connexions++;
                    }
                }
            }
            if (connexions >= 2) {
                voisins[i].push_back(j);
                voisins[j].push_back(i);
            }
        }
    }

    // Ajouter nl colonnes
    std::vector<std::set<int>> elements_Gamma(16);
    const auto& Qdata = GetData(Q);
    for (const auto& [j, pcol, val] : Qdata) {
        if (val != 0.0) {
            elements_Gamma[static_cast<int>(pcol)].insert(j);
        }
    }
    for (int p=0; p<16; ++p) {
        for (std::size_t couche=0; couche<nl; ++couche) {
            std::set<int> new_couche;
            for (auto j : elements_Gamma[p]) {
                for (auto v : voisins[j]) { 
                    new_couche.insert(v);
                }
            }
        elements_Gamma[p].insert(new_couche.begin(), new_couche.end());
        }
    }

    // Construire Gamma et R
    Mesh2DPart Gamma(16);
    CooMatrix<double> R(ne, 16);
    for (int p=0; p<16; ++p) {
        int index_local = 1;
        for (int j=0; j<ne; ++j) {
            if (elements_Gamma[p].count(j)) {
                Gamma[p].push_back(Omega[j]);
                R.push_back(j, p, index_local++);
            }
        }
    }

    // Reconstruire chaque sous-maillage avec une liste de noeuds locaux
    for (int p=0; p<16; ++p) {
        Nodes noeuds_submesh;
        for (const auto& elem : Gamma[p]) {
            for (int v=0; v<3; ++v) {
                const auto& arete = elem[v];
                if (find_index(arete, noeuds_submesh) == -1) {
                    noeuds_submesh.push_back(arete);
                }
            }
        }
        Mesh2D reconstruction(noeuds_submesh);
        for (const auto& elem : Gamma[p]) {
            reconstruction.push_back(elem);
        }
        Gamma[p] = std::move(reconstruction);
    }
    return make_pair(Gamma, R);
 }

 void
Plot(const std::vector<Mesh2D>& Sigma, const std::string& filename) {
    std::filesystem::path file = std::filesystem::path(filename);
    std::vector<std::string> tag = {"Vertices", "Edges", "Triangles", "Tetrahedra"};

    file.replace_extension(".mesh");
    std::ofstream f;
    f.open(file.c_str());

    f << "MeshVersionFormatted 3\n";
    f << "Dimension\n3\n";

    int vertices_size = 0;
    int element_size = 0;
    std::vector<Nodes> Noeuds(Sigma.size());
    for (int i=0; i<int(Sigma.size());i++) {
        Noeuds[i] = Sigma[i].nodes();
        vertices_size += Noeuds[i].size();
        element_size += Sigma[i].size();
    }

    // Vertices
    f << "Vertices\n";
    f << vertices_size << "\n";
    int k = 0;
    for (int i=0; i<int(Sigma.size()); i++) {
        for (auto& x : Noeuds[i]) {
            f << x << "\t" << i+1 << "\n";
            k++;
        }
    }
    f << "\n";

    // Elements
    f << tag[2] << "\n";
    f << element_size << "\n";

    int id, temp, offset;
    for (int i=0; i<int(Sigma.size()); i++) {
        for (const auto& elem : Sigma[i]) {
                for (std::size_t j=0; j<3;++j) {
                    id = 0;
                    temp = 0;
                    offset = 0;
                    k = 0;
                    while (id == 0 && k < int(Sigma.size())) {
                        temp = find_index(elem[j], Noeuds[k]);
                        if (temp != -1) {
                            id = temp + offset;
                        }
                        offset += Noeuds[k].size();
                        k++;
                    }
                    f << 1 + id << "\t";
                }
                f << i+1 << "\n";
        }
    }
    f << "End";
    f.close();
}
