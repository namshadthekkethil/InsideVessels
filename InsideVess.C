#include "InsideVess.h"

using namespace libMesh;
using namespace std;

Vess InsideVess::vess_i;
vector<Vess> InsideVess::vessels;
int InsideVess::vess_start;

InsideVess::InsideVess() {}

InsideVess::~InsideVess() {}

void InsideVess::read_vessel_data() {

  ifstream file_tree;
  file_tree.open("updated_vessels_Lee_data.dat");

  int vess_counter = 0;
  while (!file_tree.eof()) {
    file_tree >> vess_i.x1 >> vess_i.y1 >> vess_i.z1 >> vess_i.x2 >>
        vess_i.y2 >> vess_i.z2 >> vess_i.l >> vess_i.r1 >> vess_i.r2 >>
        vess_i.r >> vess_i.p >> vess_i.dl >> vess_i.dr;

    if (file_tree.eof())
      break;

    else {

      if (vess_i.dl == -10)
        vess_i.ter = 1;
      else if (vess_i.p == -10)
        vess_i.ter = -1;
      else
        vess_i.ter = 0;

      if (vess_i.dl != -10 && vess_i.dr != -10)
        vess_i.ter = 2;

      vessels.push_back(vess_i);

      if (vess_i.p == -10)
        vess_start = vess_counter;

      vess_counter++;
    }
  }

  for (int i = 0; i < vessels.size(); i++) {
    if (vessels[i].ter == 2) {
      vessels[vessels[i].dl].ter = 3;
      vessels[vessels[i].dr].ter = 3;
    }
  }

  file_tree.close();
}

void InsideVess::read_mesh(Mesh &mesh) {

  mesh.read("LVtrunk_heart_real_Lee_coupl.xda", NULL);
}

void InsideVess::update_inside(int n) {
  if (vessels[n].dl != -10)
    vessels[vessels[n].dl].inside = 1;
  if (vessels[n].dr != -10)
    vessels[vessels[n].dr].inside = 1;

  if (vessels[n].dl != -10)
    update_inside(vessels[n].dl);
  if (vessels[n].dr != -10)
    update_inside(vessels[n].dr);
}

void InsideVess::vessel_inside_outside(Mesh &mesh) {

  for (unsigned int i = 0; i < vessels.size(); i++) {
    vessels[i].inside = 0;
  }

  for (unsigned int i = 0; i < vessels.size(); i++) {
    cout << "i=" << i << endl;
    if (vessels[i].dl != -10 && vessels[i].dr != -10) {

      double x1 = vessels[i].x1;
      double y1 = vessels[i].y1;
      double z1 = vessels[i].z1;

      double x2 = vessels[i].x2;
      double y2 = vessels[i].y2;
      double z2 = vessels[i].z2;

      int inside_LV = 0;

      inside_elements(mesh, 0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (z1 + z2),
                      inside_LV);

      // vessels[vessels[i].dl].inside = inside_LV;
      // vessels[vessels[i].dr].inside = inside_LV;

      if (inside_LV == 1) {
        vessels[vessels[i].dl].inside = 1;
        vessels[vessels[i].dr].inside = 1;

        update_inside(vessels[i].dl);
        update_inside(vessels[i].dr);
      }
    }
  }
}

void InsideVess::inside_elements(Mesh &mesh, double xl, double yl, double zl,
                                 int &inside_LV) {
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  inside_LV = 0;

  for (; el != end_el; ++el) {
    Elem *elem = *el;
    const int elem_id = elem->id();

    double x_elem = 0.0;
    double y_elem = 0.0;
    double z_elem = 0.0;
    for (unsigned int j = 0; j < 4; j++) {

      Point pj;
      pj = elem->point(j);

      x_elem += pj(0);
      y_elem += pj(1);
      z_elem += pj(2);
    }
    x_elem /= 4.0;
    y_elem /= 4.0;
    z_elem /= 4.0;

    double h_elem = min(elem->length(0, 1), elem->length(1, 2));
    h_elem = min(h_elem, elem->length(0, 2));

    h_elem = min(h_elem, elem->length(0, 3));
    h_elem = min(h_elem, elem->length(1, 3));
    h_elem = min(h_elem, elem->length(2, 3));

    double dist_elem =
        sqrt(pow(xl - x_elem, 2) + pow(yl - y_elem, 2) + pow(zl - z_elem, 2));

    if (dist_elem < 0.5 * h_elem)
      inside_LV = 1;
  }
}

void InsideVess::writeData() {
  ofstream file_vess;
  file_vess.open("updated_vessels_Lee_inside_mod.csv", ios::out);
  file_vess << "\"x1\""
            << ",\"y1\""
            << ",\"z1\""
            << ",\"x2\""
            << ",\"y2\""
            << ",\"z2\""
            << ",\"l\""
            << ",\"r1\""
            << ",\"r2\""
            << ",\"pr\""
            << ",\"dl\""
            << ",\"dr\""
            << ",\"inside\"" << endl;
  for (int i = 0; i < vessels.size(); i++) {
    double x1_vess = vessels[i].x1;
    double y1_vess = vessels[i].y1;
    double z1_vess = vessels[i].z1;

    double x2_vess = vessels[i].x2;
    double y2_vess = vessels[i].y2;
    double z2_vess = vessels[i].z2;

    double l_vess = sqrt((x2_vess - x1_vess) * (x2_vess - x1_vess) +
                         (y2_vess - y1_vess) * (y2_vess - y1_vess) +
                         (z2_vess - z1_vess) * (z2_vess - z1_vess));
    double r1_vess = vessels[i].r; // r[vessels[i].p1];
    double r2_vess = vessels[i].r; // r[vessels[i].p2];

    file_vess << x1_vess << "," << y1_vess << "," << z1_vess << "," << x2_vess
              << "," << y2_vess << "," << z2_vess << "," << l_vess << ","
              << r1_vess << "," << r2_vess << "," << vessels[i].p << ","
              << vessels[i].dl << "," << vessels[i].dr << ","
              << vessels[i].inside << endl;
  }
  file_vess.close();

  ofstream file_vess_dat;
  file_vess_dat.open("updated_vessels_Lee_inside_data_mod.dat", ios::out);

  for (int i = 0; i < vessels.size(); i++) {
    double x1_vess = vessels[i].x1;
    double y1_vess = vessels[i].y1;
    double z1_vess = vessels[i].z1;

    double x2_vess = vessels[i].x2;
    double y2_vess = vessels[i].y2;
    double z2_vess = vessels[i].z2;

    double l_vess = sqrt((x2_vess - x1_vess) * (x2_vess - x1_vess) +
                         (y2_vess - y1_vess) * (y2_vess - y1_vess) +
                         (z2_vess - z1_vess) * (z2_vess - z1_vess));
    double r1_vess = vessels[i].r; // r[vessels[i].p1];
    double r2_vess = vessels[i].r; // r[vessels[i].p2];
    file_vess_dat << x1_vess << " " << y1_vess << " " << z1_vess << " "
                  << x2_vess << " " << y2_vess << " " << z2_vess << " "
                  << l_vess << " " << r1_vess << " " << r2_vess << " "
                  << 0.5 * (r1_vess + r2_vess) << " " << vessels[i].p << " "
                  << vessels[i].dl << " " << vessels[i].dr << " "
                  << vessels[i].inside << endl;

    // cout<<"i="<<i<<" r1="<<r1_vess<<" r2="<<r2_vess<<"
    // "<<r1_vess-r2_vess<<endl;
  }
  file_vess_dat.close();
}
