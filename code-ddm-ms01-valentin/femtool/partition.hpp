#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <filesystem>
#include <set>
#include <utility>

using Mesh2DPart = std::vector<Mesh2D>;
using FeSpace2D = FeSpace<2>;
using FeSpace2DxCoo = std::pair<FeSpace2D, CooMatrix<double>>;

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


FeSpace2DxCoo 
 Restrict(const FeSpace2D& Vh, const Mesh2D& Gamma, const std::vector<std::size_t>& tbl) {

    // Espace élément fini Uh
    FeSpace2D Uh(Gamma);

    // Maillages
    const auto& Omega = Vh.mesh();
    const auto& noeudsOmega = Omega.nodes();
    const auto& noeudsGamma = Gamma.nodes();

    // Construction de l'indexage global -> local
    std::vector<int> Gamma2Omega(noeudsGamma.size(), -1);
    for (std::size_t k=0; k<Gamma.size(); ++k) {
        auto& elemGamma = Gamma[k];
        auto& elemOmega = Omega[tbl[k]];
        for (int v=0; v<3; ++v) {
            int idGamma = find_index(elemGamma[v], noeudsGamma); // indice du v-ème sommet de Gamma
            int idOmega = find_index(elemOmega[v], noeudsOmega); // indice du v-ème sommet correspondant
            Gamma2Omega[idGamma] = idOmega;
        }        
    }

    // Construction de P
    std::size_t nUh = dim(Uh);
    std::size_t nVh = dim(Vh);
    CooMatrix<double> P(nUh, nVh);
    // Éléments de Uh
    for (std::size_t k=0; k<Uh.size(); ++k) {
        auto& cell = Uh[k];
        for (std::size_t loc=0; loc<3; ++loc) {
            std::size_t dofU = cell[loc]; // dof global dans Uh
            std::size_t noeudGamma = loc; // indice du noeud local dans Gamma
            int noeud0 = Gamma2Omega[noeudGamma]; // on retrouve quel noeud d'Omega correspond à ce noeud de Omega
            if (noeud0 < 0) {
                continue; // si pas de correspondance on peut sauter ce noeud (Gamma2Omega initialisée à -1)
            }

            // Éléments de Vh: on cherche dans Vh quel DoF correspond à noeud0
            for (std::size_t l=0; l<Vh.size(); ++l) {
                auto& cell2 = Vh[l];
                for (std::size_t loc2=0; loc2<3; ++loc2) {
                    if (Close(Omega[l][loc2], noeudsOmega[noeud0])) {
                        std::size_t dofV = cell2[loc2]; // numéro de dof dans Vh pour ce noeud
                        P.push_back(dofU, dofV, 1.0);
                    }
                }
            }
        }
    }
    return std::make_pair(Uh, P);
}

std::vector<FeSpace2DxCoo> 
 Partition4(const FeSpace2D& Vh, const std::size_t& nl) {

    Mesh2D Omega = Vh.mesh();
    auto [Sigma, Q] = Partition4(Omega, nl);

    std::vector<std::vector<std::size_t>> tbl_global(Sigma.size());
    for (std::size_t p=0; p<Sigma.size(); ++p) {
        tbl_global[p].resize(Sigma[p].size());
        for (std::size_t j=0; j<Sigma[p].size(); ++j) {
            auto& elemSigma = Sigma[p][j];
            for (std::size_t k=0; k<Omega.size(); ++k) {
                if (elemSigma == Omega[k]) {
                    tbl_global[p][j] = k;
                    break;
                }
            }
        }
    }

    std::vector<FeSpace2DxCoo> UP(Sigma.size());
    for (std::size_t p=0; p<Sigma.size(); ++p) {
        UP[p] = Restrict(Vh, Sigma[p], tbl_global[p]);
    }
    return UP;
}


std::vector<FeSpace2DxCoo> 
 Partition16(const FeSpace2D& Vh, const std::size_t& nl) {

    Mesh2D Omega = Vh.mesh();
    auto [Sigma, Q] = Partition4(Omega, nl);

    std::vector<std::vector<std::size_t>> tbl_global(Sigma.size());
    for (std::size_t p=0; p<Sigma.size(); ++p) {
        tbl_global[p].resize(Sigma[p].size());
        for (std::size_t j=0; j<Sigma[p].size(); ++j) {
            auto& elemSigma = Sigma[p][j];
            for (std::size_t k=0; k<Omega.size(); ++k) {
                if (elemSigma == Omega[k]) {
                    tbl_global[p][j] = k;
                    break;
                }
            }
        }
    }

    std::vector<FeSpace2DxCoo> UP(Sigma.size());
    for (std::size_t p=0; p<Sigma.size(); ++p) {
        UP[p] = Restrict(Vh, Sigma[p], tbl_global[p]);
    }
    return UP;
}

