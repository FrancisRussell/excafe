#ifndef SIMPLE_CFD_FORMS_FORM_EVALUATOR_HPP
#define SIMPLE_CFD_FORMS_FORM_EVALUATOR_HPP

#include <cstddef>
#include <simple_cfd/numeric/tensor.hpp>
#include <simple_cfd/cell_vertices.hpp>
#include <simple_cfd/vertex.hpp>
#include <simple_cfd/general_cell.hpp>
#include <simple_cfd/cell_manager.hpp>
#include <simple_cfd/mesh_entity.hpp>
#include <simple_cfd/capture/capture_fwd.hpp>
#include <simple_cfd/capture/evaluation/evaluation_fwd.hpp>
#include "linear_form.hpp"
#include "form_evaluator_visitor.hpp"

namespace cfd
{

namespace forms
{

template<std::size_t D>
class FormEvaluator
{
private:
  static const std::size_t dimension = D;
  typedef typename CellManager::ref<dimension>::general cell_ref_t;
  cell_ref_t cell;
  const Scenario<dimension>* scenario;
  const detail::ExpressionValues<dimension>* values;
  LinearForm form;

public:
  FormEvaluator(const cell_ref_t _cell, const Scenario<dimension>& _scenario, 
    const detail::ExpressionValues<dimension>& _values, const LinearForm& _form) :
    cell(_cell), scenario(&_scenario), values(&_values), form(_form)
  {
  }

  FormEvaluator(const FormEvaluator& f) : cell(f.cell), scenario(f.scenario),
    values(f.values), form(f.form)
  {
  }

  Tensor<dimension> evaluate(const CellVertices<dimension>& vertices, 
    const MeshEntity& localEntity, const vertex<dimension>& v, const Dof<dimension>& dof) const
  {
    FormEvaluatorVisitor<dimension> visitor(cell, *scenario, *values, vertices, localEntity, v, dof);
    form.accept(visitor);
    return visitor.getResult();
  }
};

}

}

#endif
