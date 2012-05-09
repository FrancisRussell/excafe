#include <excafe/cse/kcm.hpp>
#include <excafe/cse/sop.hpp>
#include <excafe/cse/cube.hpp>
#include <excafe/cse/biclique.hpp>
#include <excafe/cse/biclique_search.hpp>
#include <excafe/cse/polynomial_index.hpp>
#include <excafe/cse/new_literal_creator.hpp>
#include <excafe/exception.hpp>
#include <map>
#include <utility>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>

namespace excafe
{

namespace cse
{

void KCM::addPolynomial(const PolynomialIndex& polynomialID)
{
  const SOP& sop = sops[polynomialID];
  const SOP::cokernel_set_t coKernels = sop.getCoKernels();
  BOOST_FOREACH(const Cube& coKernel, coKernels)
  {
    const vertex_descriptor coKernelVertex = addCoKernel(polynomialID, coKernel);
    const SOP kernel = sop/coKernel;

    for(SOP::const_iterator iter = kernel.begin(); iter != kernel.end(); ++iter)
    {
      const vertex_descriptor cubeVertex = addCube(*iter);
      const std::size_t termID = kernel.getTermNumber(iter);
      const std::pair<edge_descriptor, bool> edgePair = add_edge(coKernelVertex, cubeVertex, graph);

      if (!edgePair.second)
        CFD_EXCEPTION("Attemped to insert duplicate edge into KCM. This should never happen.");

      const edge_descriptor edge = edgePair.first;
      put(term_id(), graph, edge, termID);
    }
  }
}

KCM::KCM(NewLiteralCreator& _literalCreator) : 
literalCreator(_literalCreator), sops(literalCreator.getSOPMap())
{
  BOOST_FOREACH(const SOPMap::value_type mapping, sops)
  {
    addPolynomial(mapping.first);
  }
}

KCM::vertex_descriptor KCM::addCoKernel(const PolynomialIndex& polynomialID, const Cube& coKernel)
{
  const vertex_descriptor v = add_vertex(graph);
  put(is_cube(), graph, v, false);
  put(polynomial_id(), graph, v, polynomialID);
  put(term_cube(), graph, v, coKernel);
  put(mul_count(), graph, v, coKernel.numMultiplies(literalCreator));
  put(is_unit(), graph, v, coKernel.isUnit(literalCreator));
  put(is_numeric(), graph, v, coKernel.isNumeric(literalCreator));
  put(has_coefficient(), graph, v, coKernel.hasCoefficient(literalCreator));
  return v;
}

KCM::vertex_descriptor KCM::addCube(const Cube& c)
{
  const std::map<Cube, vertex_descriptor>::const_iterator vertexIter = cubeVertices.find(c);

  if (vertexIter != cubeVertices.end())
  {
    return vertexIter->second;
  }
  else
  {
    const vertex_descriptor v = add_vertex(graph);
    cubeVertices.insert(std::make_pair(c, v));
    put(is_cube(), graph, v, true);
    put(term_cube(), graph, v, c);
    put(mul_count(), graph, v, c.numMultiplies(literalCreator));
    put(is_unit(), graph, v, c.isUnit(literalCreator));
    put(is_numeric(), graph, v, c.isNumeric(literalCreator));
    put(has_coefficient(), graph, v, c.hasCoefficient(literalCreator));
    return v;
  }
}

void KCM::orderCubes()
{
  unsigned id=0;
  BOOST_FOREACH(const vertex_descriptor& vertex, vertices(graph))
  {
    if (get(is_cube(), graph, vertex))
    {
      put(cube_ordering(), graph, vertex, std::make_pair(out_degree(vertex, graph), id));
      ++id;
    }
  }
}

PolynomialIndex KCM::addPolynomial(const SOP& sop)
{
  const PolynomialIndex index = sops.addSOP(sop);
  addPolynomial(index);
  return index;
}

bool KCM::factorise()
{
  orderCubes();

  std::priority_queue<biclique_search_t, std::vector<biclique_search_t>, BicliqueSearchComparator> queue;
  addSearchSpaces(queue);
  biclique_t best(graph);

  while(!queue.empty())
  {
    const biclique_search_t bs = queue.top();
    queue.pop();

    if (bs.getMaximalValue() < best.getValue())
    {
      break;
    }
    else
    {
      //bs.print();
      //std::cout << std::endl;

      if (!bs.isFinished())
      {
        const std::pair<biclique_search_t, biclique_search_t> pair = bs.split();

        if (best.getValue() <= pair.first.getMaximalValue())
          queue.push(pair.first);

        if (best.getValue() <= pair.second.getMaximalValue())
          queue.push(pair.second);
      }

      const bool newBest = best.getValue() < bs.getValue()
                           || (best.getValue() == bs.getValue()
                               && best.numCubes()*best.numCoKernels() > bs.numCubes()*bs.numCoKernels());
      if (newBest)
      {
        best = bs.getBiclique();
        std::cout << "New best score: " << best.getValue() << std::endl;
      }
    }
  }

  if (best.getValue() > 0)
  {
    removeBiclique(best);
    return true;
  }
  else
  {
    return false;
  }
}

std::size_t KCM::numEdges() const
{
  return num_edges(graph);
}

std::size_t KCM::numCubes() const
{
  std::size_t count = 0;
  BOOST_FOREACH(const vertex_descriptor& v, vertices(graph))
  {
    if (get(is_cube(), graph, v))
      ++count;
  }
  return count;
}

std::size_t KCM::numCoKernels() const
{
  return num_vertices(graph) - numCubes();
}

std::size_t KCM::numAdditions() const
{
  std::size_t result=0;
  BOOST_FOREACH(const SOPMap::value_type& mapping, sops)
    result += mapping.second.numAdditions(literalCreator);

  return result;
}

std::size_t KCM::numMultiplies() const
{
  std::size_t result=0;
  BOOST_FOREACH(const SOPMap::value_type& mapping, sops)
    result += mapping.second.numMultiplies(literalCreator);

  return result;
}

void KCM::updateGraph(const Biclique<graph_t>& biclique)
{
  const std::set<PolynomialIndex> modifiedPolynomials(biclique.getModifiedPolynomials());

  typedef boost::graph_traits<graph_t>::vertex_iterator vertex_iter;
  vertex_iter vi, viEnd;

  boost::tie(vi, viEnd) = vertices(graph);
  while(vi != viEnd)
  {
    const vertex_iter viNext = boost::next(vi);
    if (!get(is_cube(), graph, *vi))
    {
      if (modifiedPolynomials.find(get(polynomial_id(), graph, *vi)) != modifiedPolynomials.end())
      {
        clear_vertex(*vi, graph);
        remove_vertex(*vi, graph);
      }
    }

    vi = viNext;
  }

  BOOST_FOREACH(const PolynomialIndex& polynomialID, modifiedPolynomials)
  {
    addPolynomial(polynomialID);
  }
}

void KCM::removeBiclique(const Biclique<graph_t>& biclique)
{
  const SOP newSOP = biclique.getSOP();
  const PolynomialIndex newSOPIndex = addPolynomial(newSOP);
  const unsigned literal = literalCreator.getLiteralID(newSOPIndex);

  // Rewrite polynomials
  biclique.rewritePolynomials(addCube(literal), sops);

  // Update graph
  updateGraph(biclique);
}

}

}
