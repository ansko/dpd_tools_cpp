#ifndef OPTIONS_HPP
#define OPTIONS_HPP


#include <iostream>
#include <fstream>
#include <map>
#include <string>


class Options
{
public:
    size_t md_soft_atoms;
    size_t platelet_radius;
    size_t tail_length;
    size_t stacking;
    size_t polymerization;
    size_t mmt_atom_type;
    size_t mmt_edge_bond_type;
    size_t mmt_diagonal_bond_type;
    size_t modifier_head_atom_type;
    size_t modifier_tail_atom_type;
    size_t head_tail_type;
    size_t tail_tail_type;
    size_t polymer_atom_type;
    size_t polymer_bond_type;
    size_t platelet_edge;
    size_t modifiers_count_preset = 0;
    float lx;
    float ly;
    float lz;
    float mmt_real_thickness;
    float planar_expansion_coeff;
    float bead_charge;
    float real_interlayer;
    float dpd_rho;
    float lj_bead_radius;
    float too_close_threshold_mmt;
    float too_close_threshold_soft;
    float modifier_head_tail_bond_length;
    float modifier_tail_tail_bond_length;
    float polymer_bond_length;

    Options(std::string options_fname)
      {
        std::string option_name;
        std::ifstream ofs(options_fname);
        while (ofs >> option_name)
          {
            if (option_name == "lx")
                ofs >> this->lx;
            else if (option_name == "ly")
                ofs >> this->ly;
            else if (option_name == "lz")
                ofs >> this->lz;
            else if (option_name == "mmt_real_thickness")
                ofs >> this->mmt_real_thickness;
            else if (option_name == "planar_expansion_coeff")
                ofs >> this->planar_expansion_coeff;
            else if (option_name == "bead_charge")
                ofs >> this->bead_charge;
            else if (option_name == "real_interlayer")
                ofs >> this->real_interlayer;
            else if (option_name == "dpd_rho")
                ofs >> this->dpd_rho;
            else if (option_name == "md_soft_atoms")
                ofs >> this->md_soft_atoms;
            else if (option_name == "platelet_radius")
                ofs >> this->platelet_radius;
            else if (option_name == "tail_length")
                ofs >> this->tail_length;
            else if (option_name == "stacking")
                ofs >> this->stacking;
            else if (option_name == "polymerization")
                ofs >> this->polymerization;
            else if (option_name == "lj_bead_radius")
                ofs >> this->lj_bead_radius;
            else if (option_name == "mmt_atom_type")
                ofs >> this->mmt_atom_type;
            else if (option_name == "mmt_edge_bond_type")
                ofs >> this->mmt_edge_bond_type;
            else if (option_name == "mmt_diagonal_bond_type")
                ofs >> this->mmt_diagonal_bond_type;
            else if (option_name == "modifier_head_atom_type")
                ofs >> this->modifier_head_atom_type;
            else if (option_name == "modifier_tail_atom_type")
                ofs >> this->modifier_tail_atom_type;
            else if (option_name == "head_tail_type")
                ofs >> this->head_tail_type;
            else if (option_name == "tail_tail_type")
                ofs >> this->tail_tail_type;
            else if (option_name == "too_close_threshold_mmt")
                ofs >> this->too_close_threshold_mmt;
            else if (option_name == "too_close_threshold_soft")
                ofs >> this->too_close_threshold_soft;
            else if (option_name == "modifier_head_tail_bond_length")
                ofs >> this->modifier_head_tail_bond_length;
            else if (option_name == "modifier_tail_tail_bond_length")
                ofs >> this->modifier_tail_tail_bond_length;
            else if (option_name == "polymer_bond_length")
                ofs >> this->polymer_bond_length;
            else if (option_name == "polymer_atom_type")
                ofs >> this->polymer_atom_type;
            else if (option_name == "polymer_bond_type")
                ofs >> this->polymer_bond_type;
            else if (option_name == "platelet_edge")
                ofs >> this->platelet_edge;
            else if (option_name == "modifiers_count_preset")
                ofs >> this->modifiers_count_preset;
          }
      };

    void print_options(bool verbose=false)
      {
        std::cout << "lx = " << this->lx;
        if (verbose)
            std::cout << " Cell length along x axis\n";
        else
            std::cout << std::endl;
        std::cout << "ly = " << this->ly;
        if (verbose)
            std::cout << " Cell length along y axis\n";
        else
            std::cout << std::endl;
        std::cout << "lz = " << this->lz;
        if (verbose)
            std::cout << " Cell length along z axis\n";
        else
            std::cout << std::endl;
        std::cout << "mmt_real_thickness = " << this->mmt_real_thickness;
        if (verbose)
            std::cout << " Thickness of a single MMT platelet (in Angstroms)\n";
        else
            std::cout << std::endl;
        std::cout << "planar_expansion_coeff = " << this->planar_expansion_coeff;
        if (verbose)
            std::cout << " Cell increase coeff in xy plane around MMT platelet\n";
        else
            std::cout << std::endl;
        std::cout << "bead_charge = " << this->bead_charge;
        if (verbose)
            std::cout << " Charge of every charged bead (in e)\n";
        else
            std::cout << std::endl;
        std::cout << "real_interlayer = " << this->real_interlayer;
        if (verbose)
            std::cout << " Interlayer thickness in MMT (in Angstroms)\n";
        else
            std::cout << std::endl;
        std::cout << "dpd_rho = " << this->dpd_rho;
        if (verbose)
            std::cout << " DPD density (beads count / volume)\n";
        else
            std::cout << std::endl;
        std::cout << "md_soft_atoms = " << this->md_soft_atoms;
        if (verbose)
            std::cout << " Count of soft atoms in composite from our MD\n";
        else
            std::cout << std::endl;
        std::cout << "platelet_radius = " << this->platelet_radius;
        if (verbose)
            std::cout << " Radius of MMT platelet (in beads)\n";
        else
            std::cout << std::endl;
        std::cout << "tail_length = " << this->tail_length;
        if (verbose)
            std::cout << " Length of modifier's tail (in beads)\n";
        else
            std::cout << std::endl;
        std::cout << "stacking = " << this->stacking;
        if (verbose)
            std::cout << " MMT platelets count in the stack\n";
        else
            std::cout << std::endl;
        std::cout << "polymerization = " << this->polymerization;
        if (verbose)
            std::cout << " Polymerization degree (in monomers)\n";
        else
            std::cout << std::endl;
        std::cout << "lj_bead_radius = " << this->lj_bead_radius;
        if (verbose)
            std::cout << " Bead radius in lj units\n";
        else
            std::cout << std::endl;
        std::cout << "mmt_atom_type = " << this->mmt_atom_type;
        if (verbose)
            std::cout << " Type of all MMT atoms\n";
        else
            std::cout << std::endl;
        std::cout << "mmt_edge_bond_type = " << this->mmt_edge_bond_type;
        if (verbose)
            std::cout << " Type of MMT edge bonds\n";
        else
            std::cout << std::endl;
        std::cout << "mmt_diagonal_bond_type = " << this->mmt_diagonal_bond_type;
        if (verbose)
            std::cout << " Type of MMT diagonal bonds\n";
        else
            std::cout << std::endl;
        std::cout << "modifier_head_atom_type = " << this->modifier_head_atom_type;
        if (verbose)
            std::cout << " Type of modifier head atom\n";
        else
            std::cout << std::endl;
        std::cout << "modifier_tail_atom_type = " << this->modifier_tail_atom_type;
        if (verbose)
            std::cout << " Type of modifier tail atom\n";
        else
            std::cout << std::endl;
        std::cout << "head_tail_type = " << this->head_tail_type;
        if (verbose)
            std::cout << " Type of modifier head-tail bond\n";
        else
            std::cout << std::endl;
        std::cout << "tail_tail_type = " << this->tail_tail_type;
        if (verbose)
            std::cout << " Type of modifier tail-tail bond\n";
        else
            std::cout << std::endl;
        std::cout << "too_close_threshold_mmt = " << this->too_close_threshold_mmt;
        if (verbose)
            std::cout << " Threshold when beads are too close to mmt\n";
        else
            std::cout << std::endl;
        std::cout << "too_close_threshold_soft = "
            << this->too_close_threshold_soft;
        if (verbose)
            std::cout << " Threshold when beads are too close to soft phase\n";
        else
            std::cout << std::endl;
        std::cout << "modifier_head_tail_bond_length = "
            << this->modifier_head_tail_bond_length;
        if (verbose)
            std::cout << " Length of the bond between head and tail in modifier\n";
        else
            std::cout << std::endl;
        std::cout << "modifier_tail_tail_bond_length = "
            << this->modifier_tail_tail_bond_length;
        if (verbose)
            std::cout << " Length of the bond between tails in modifier\n";
        else
            std::cout << std::endl;
        std::cout << "polymer_bond_length = " << this->polymer_bond_length;
        if (verbose)
            std::cout << " Length of the bond in polymer\n";
        else
            std::cout << std::endl;
        std::cout << "polymer_atom_type = "  << this->polymer_atom_type;
        if (verbose)
            std::cout << " Type of all atoms in polymer\n";
        else
            std::cout << std::endl;
        std::cout << "polymer_bond_type = "  << this->polymer_bond_type;
        if (verbose)
            std::cout << " Type of all bonds in polymer\n";
        else
            std::cout << std::endl;
        std::cout << "platelet_edge = "  << this->platelet_edge;
        if (verbose)
            std::cout << " Size of periodic MMT platelet\n";
        else
            std::cout << std::endl;
        std::cout << "modifiers_count_preset = " << this->modifiers_count_preset;
        if (verbose)
            std::cout << " Preset number of modifiers\n";
        else
            std::cout << std::endl;
      }
};


#endif  // OPTIONS_HPP
