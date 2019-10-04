#ifndef PARAMETERS_CLASSES_HPP
#define PARAMETERS_CLASSES_HPP


struct AddMmtCircularParameters
{
public:
    AddMmtCircularParameters(OptionsParser &o, float x, float y, float z,
        size_t charged_count)
    : x(x), y(y), z(z), charged_count(charged_count),
      lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay")),
      bead_charge(o.get<float>("bead_charge")),
      platelet_closing(o.get<float>("platelet_closing")),
      platelet_radius(o.get<size_t>("platelet_radius")),
      mmt_atom_type(o.get<size_t>("mmt_atom_type")),
      mmt_edge_bond_type(o.get<size_t>("mmt_edge_bond_type")),
      mmt_diagonal_bond_type(o.get<size_t>("mmt_diagonal_bond_type"))
    {
    }

    float x;               // x of the circluar platelet center
    float y;               // y of the circluar platelet center
    float z;               // z of the circluar platelet center
    size_t charged_count;  // number of platelet atoms carrying non zero charge

    float lj_bead_radius_clay;
    float bead_charge;
    float platelet_closing;
    size_t platelet_radius;
    size_t mmt_atom_type;
    size_t mmt_edge_bond_type;
    size_t mmt_diagonal_bond_type;
};


struct AddMmtPeriodicParameters
{
public:
    AddMmtPeriodicParameters(OptionsParser &o, float z=0, size_t
        charged_count=0, std::string mode="isolated")
    : z(z), charged_count(charged_count),
      lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay")),
      bead_charge(o.get<float>("bead_charge")),
      platelet_edge(o.get<size_t>("platelet_edge")),
      mmt_atom_type(o.get<size_t>("mmt_atom_type")),
      mmt_edge_bond_type(o.get<size_t>("mmt_edge_bond_type")),
      mmt_diagonal_bond_type(o.get<size_t>("mmt_diagonal_bond_type")),
      platelet_closing(o.get<float>("platelet_closing"))
    {}

    float z;               // z of the circluar platelet center
    size_t charged_count;  // number of platelet atoms carrying non zero charge
    float lj_bead_radius_clay;
    float platelet_closing;
    float bead_charge;
    size_t platelet_edge;
    size_t mmt_atom_type;
    size_t mmt_edge_bond_type;
    size_t mmt_diagonal_bond_type;
    std::string mode;
};


struct AddModifierGalleryParameters
{
public:
    AddModifierGalleryParameters (OptionsParser &o, float top, float bottom,
                                  std::string mode)
    : top(top), bottom(bottom), mode(mode),
      modifier_head_tail_bond_length(
          o.get<float>("modifier_head_tail_bond_length")),
      modifier_tail_tail_bond_length(
          o.get<float>("modifier_tail_tail_bond_length")),
      lj_bead_radius_soft(o.get<float>("lj_bead_radius_soft")),
      lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay")),
      too_close_threshold_mmt(o.get<float>("too_close_threshold_mmt")),
      too_close_threshold_soft(o.get<float>("too_close_threshold_soft")),
      bead_charge(o.get<float>("bead_charge")),
      platelet_closing(o.get<float>("platelet_closing")),
      tail_length(o.get<size_t>("tail_length")),
      modifier_head_atom_type(o.get<size_t>("modifier_head_atom_type")),
      modifier_tail_atom_type(o.get<size_t>("modifier_tail_atom_type")),
      head_tail_type(o.get<size_t>("head_tail_type")),
      tail_tail_type(o.get<size_t>("tail_tail_type"))//,
    {
      if (mode == "isolated")
        {
          this->platelet_radius = o.get<size_t>("platelet_radius");
        }
      else if (mode == "periodic")
        ;
      else
        {
          throw std::exception();
        }
    }

    float top;     // top z-coordintate of the gallery
    float bottom;  // bottom z-coordinate of the gallery
    float modifier_head_tail_bond_length;
    float modifier_tail_tail_bond_length;
    float lj_bead_radius_soft;
    float lj_bead_radius_clay;
    float too_close_threshold_mmt;
    float too_close_threshold_soft;
    float bead_charge;
    float platelet_closing;
    size_t tail_length;
    size_t modifier_head_atom_type;
    size_t modifier_tail_atom_type;
    size_t head_tail_type;
    size_t tail_tail_type;
    size_t platelet_radius;
    std::string mode;
};


struct AddPolymerParallelParameters
{
public:
    AddPolymerParallelParameters(OptionsParser &o)
    : polymer_bond_length(o.get<float>("polymer_bond_length")),
      lj_bead_radius_soft(o.get<float>("lj_bead_radius_soft")),
      too_close_threshold_mmt(o.get<float>("too_close_threshold_mmt")),
      too_close_threshold_soft(o.get<float>("too_close_threshold_soft")),
      polymer_atom_type(o.get<size_t>("polymer_atom_type")),
      polymer_bond_type(o.get<size_t>("polymer_bond_type")),
      polymerization(o.get<size_t>("polymerization"))
    {}

    float polymer_bond_length;
    float lj_bead_radius_soft;
    float too_close_threshold_mmt;
    float too_close_threshold_soft;
    size_t polymer_atom_type;
    size_t polymer_bond_type;
    size_t polymerization;
};


#endif  // PARAMETERS_CLASSES_HPP include guard
