// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "abstract_boundary_operator.hpp"
#include "discrete_boundary_operator.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <boost/type_traits/is_complex.hpp>
#include <stdexcept>
#include <tbb/atomic.h>

namespace Bempp
{

namespace
{
    // Used to generate unique labels of anonymous operators. Automatically
    // set to zero at program startup.
    // Cannot be declared static in AbstractBoundaryOperator, since it must be
    // shared across different template instantiations.
    tbb::atomic<int> s_anonymousOperatorCounter;
}

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperator<BasisFunctionType, ResultType>::
AbstractBoundaryOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                         const shared_ptr<const Space<BasisFunctionType> >& range,
                         const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                         const std::string& label,
                         int symmetry) :
    m_domain(domain), m_range(range), m_dualToRange(dualToRange),
    m_label(label), m_symmetry(symmetry)
{
    if (!m_domain)
        throw std::invalid_argument(
                "AbstractBoundaryOperator::AbstractBoundaryOperator(): "
                "domain must not be null");
    if (!m_range)
        throw std::invalid_argument(
                "AbstractBoundaryOperator::AbstractBoundaryOperator(): "
                "range must not be null");
    if (!m_dualToRange)
        throw std::invalid_argument(
                "AbstractBoundaryOperator::AbstractBoundaryOperator(): "
                "dualToRange must not be null");
    if (&m_range->grid() != &m_dualToRange->grid())
        throw std::invalid_argument(
                "AbstractBoundaryOperator::AbstractBoundaryOperator(): "
                "range and dualToRange must be defined on the same grid");

    if (m_label.empty()) {
        int i = ++s_anonymousOperatorCounter;
        m_label = "Op" + toString(i);
    }

    if (m_symmetry & AUTO_SYMMETRY)
        m_symmetry = NO_SYMMETRY;
    // For real operators Hermitian and symmetric are equivalent
    if (!boost::is_complex<ResultType>() &&
            (m_symmetry & (SYMMETRIC | HERMITIAN)))
        m_symmetry |= SYMMETRIC | HERMITIAN;
}

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperator<BasisFunctionType, ResultType>::~AbstractBoundaryOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const AbstractBoundaryOperatorId>
AbstractBoundaryOperator<BasisFunctionType, ResultType>::id() const
{
    return shared_ptr<const AbstractBoundaryOperatorId>();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
AbstractBoundaryOperator<BasisFunctionType, ResultType>::assembleWeakForm(
        const Context<BasisFunctionType, ResultType>& context) const
{
    if (!m_domain->dofsAssigned())
        throw std::runtime_error(
                "AbstractBoundaryOperator::assembleWeakForm(): "
                "assignDofs() has not been called on the domain space");
    if (!m_dualToRange->dofsAssigned())
        throw std::runtime_error(
                "AbstractBoundaryOperator::assembleWeakForm(): "
                "assignDofs() has not been called on the space dual to range");

    return this->assembleWeakFormImpl(context);
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType> >
AbstractBoundaryOperator<BasisFunctionType, ResultType>::domain() const
{
    return m_domain;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType> >
AbstractBoundaryOperator<BasisFunctionType, ResultType>::range() const
{
    return m_range;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType> >
AbstractBoundaryOperator<BasisFunctionType, ResultType>::dualToRange() const
{
    return m_dualToRange;
}

template <typename BasisFunctionType, typename ResultType>
std::string
AbstractBoundaryOperator<BasisFunctionType, ResultType>::label() const
{
    return m_label;
}

template <typename BasisFunctionType, typename ResultType>
int
AbstractBoundaryOperator<BasisFunctionType, ResultType>::symmetry() const
{
    return m_symmetry;
}

template <typename BasisFunctionType, typename ResultType>
void
AbstractBoundaryOperator<BasisFunctionType, ResultType>::collectDataForAssemblerConstruction(
        const AssemblyOptions& options,
        shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
        shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
        shared_ptr<GeometryFactory>& testGeometryFactory,
        shared_ptr<GeometryFactory>& trialGeometryFactory,
        shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
        shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
        shared_ptr<Fiber::OpenClHandler>& openClHandler,
        bool& cacheSingularIntegrals) const
{
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;
    typedef LocalAssemblerConstructionHelper Helper;

    // Collect grid data
    Helper::collectGridData(m_dualToRange->grid(),
                            testRawGeometry, testGeometryFactory);
    if (&m_dualToRange->grid() == &m_domain->grid()) {
        trialRawGeometry = testRawGeometry;
        trialGeometryFactory = testGeometryFactory;
    } else
        Helper::collectGridData(m_domain->grid(),
                                trialRawGeometry, trialGeometryFactory);

    // Construct the OpenClHandler
    Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                              testRawGeometry, trialRawGeometry, openClHandler);

    // Get pointers to test and trial bases of each element
    Helper::collectBases(*m_dualToRange, testBases);
    if (m_dualToRange == m_domain)
        trialBases = testBases;
    else
        Helper::collectBases(*m_domain, trialBases);

    cacheSingularIntegrals = options.isSingularIntegralCachingEnabled();
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);

} // namespace Bempp