

#ifndef bempp_general_elementary_singular_integral_operator_hpp
#define bempp_general_elementary_singular_integral_operator_hpp

#include "elementary_singular_integral_operator.hpp"

//Peter:
#include "abstract_boundary_operator.hpp"
#include "dense_global_assembler.hpp"

namespace Bempp
{
template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
class GeneralElementarySingularIntegralOperator :
        public ElementarySingularIntegralOperator<
            BasisFunctionType_, KernelType_, ResultType_>
{
    typedef ElementarySingularIntegralOperator<
    BasisFunctionType_, KernelType_, ResultType_> Base;
public:
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \brief Type of the values of kernel functions. */
    typedef typename Base::KernelType KernelType;
    /** \copydoc ElementaryIntegralOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ElementaryIntegralOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ElementaryIntegralOperator::CollectionOfShapesetTransformations */
    typedef typename Base::CollectionOfShapesetTransformations
    CollectionOfShapesetTransformations;
    /** \copydoc ElementaryIntegralOperator::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc ElementaryIntegralOperator::CollectionOfKernels */
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    /** \copydoc ElementaryIntegralOperator::TestKernelTrialIntegral */
    typedef typename Base::TestKernelTrialIntegral TestKernelTrialIntegral;

    template <typename KernelFunctor, typename TestTransformationsFunctor, typename TrialTransformationsFunctor, typename IntegrandFunctor>
    GeneralElementarySingularIntegralOperator(
            const shared_ptr<const Space<BasisFunctionType_> >& domain,
            const shared_ptr<const Space<BasisFunctionType_> >& range,
            const shared_ptr<const Space<BasisFunctionType_> >& dualToRange,
            const std::string& label,
            int symmetry,
            const KernelFunctor& kernelFunctor,
            const TestTransformationsFunctor& testTransformationsFunctor,
            const TrialTransformationsFunctor& trialTransformationsFunctor,
            const IntegrandFunctor& integrandFunctor);

    /** \overload
     *
     *  This constructor takes the same arguments as the preceding one
     *  except for the \p integral argument, which should be a shared pointer
     *  to an instance of (a subclass of) Fiber::TestKernelTrialIntegral.
     */
    template <typename KernelFunctor,typename TestTransformationsFunctor,typename TrialTransformationsFunctor>
    GeneralElementarySingularIntegralOperator(
            const shared_ptr<const Space<BasisFunctionType_> >& domain,
            const shared_ptr<const Space<BasisFunctionType_> >& range,
            const shared_ptr<const Space<BasisFunctionType_> >& dualToRange,
            const std::string& label,
            int symmetry,
            const KernelFunctor& kernelFunctor,
            const TestTransformationsFunctor& testTransformationsFunctor,
            const TrialTransformationsFunctor& trialTransformationsFunctor,
            const shared_ptr<Fiber::TestKernelTrialIntegral<
            BasisFunctionType_, KernelType_, ResultType_> >& integral);


    /** \overload
     *
     *  This constructor takes the same first five arguments as the preceding
     *  ones, but the last four arguments should be shared pointers to
     *  instances of Fiber::CollectionOfKernels,
     *  Fiber::CollectionOfShapesetTransformations and
     *  Fiber::TestKernelTrialIntegral.
     */
    GeneralElementarySingularIntegralOperator(
        const shared_ptr<const Space<BasisFunctionType_> >& domain,
        const shared_ptr<const Space<BasisFunctionType_> >& range,
        const shared_ptr<const Space<BasisFunctionType_> >& dualToRange,
        const std::string& label,
        int symmetry,
        const shared_ptr<Fiber::CollectionOfKernels<KernelType_> >& kernels,
        const shared_ptr<Fiber::CollectionOfShapesetTransformations<CoordinateType> >&
        testTransformations,
        const shared_ptr<Fiber::CollectionOfShapesetTransformations<CoordinateType> >&
        trialTransformations,
        const shared_ptr<Fiber::TestKernelTrialIntegral<
        BasisFunctionType_, KernelType_, ResultType_> >& integral);

    virtual const CollectionOfKernels& kernels() const
    { return *m_kernels; }
    virtual const CollectionOfShapesetTransformations& testTransformations() const
    { return *m_testTransformations; }
    virtual const CollectionOfShapesetTransformations& trialTransformations() const
    { return *m_trialTransformations; }
    virtual const TestKernelTrialIntegral& integral() const
    { return *m_integral; }


//Peter:
//shared_ptr<const DiscreteBoundaryOperator<ResultType> > weakFormPeter(std::string str, const Context<BasisFunctionType, ResultType> context, std::unique_ptr<arma::Col<ResultType> > solV, std::unique_ptr<std::vector<ResultType> > rhsV, std::unique_ptr<arma::Mat<ResultType> > wm) const
shared_ptr<const DiscreteBoundaryOperator<ResultType> > weakFormPeter(std::string str, const Context<BasisFunctionType, ResultType> context, arma::Col<ResultType> * solV, std::vector<ResultType> * rhsV, arma::Mat<ResultType> * wm) const
{
//	std::cout << "saoidjfid" << std::endl;
	typedef AbstractBoundaryOperator<BasisFunctionType, ResultType> Base;
	typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssembler;
//std::auto_ptr<LocalAssembler> assembler = this->makeAssembler(*context.quadStrategy(), context.assemblyOptions());
	std::unique_ptr<LocalAssembler> assembler = this->makeAssembler(*context.quadStrategy(), context.assemblyOptions());
//	std::cout << context.assemblyOptions().assemblyMode() <<  " iosujhdfoiakd " << AssemblyOptions::DENSE << " asioufhd " << AssemblyOptions::ACA << std::endl;
	switch (context.assemblyOptions().assemblyMode()) {
	case AssemblyOptions::DENSE: {
//		std::cout << "oiasujhdfid" << std::endl;
//        return shared_ptr<DiscreteBoundaryOperator<ResultType> >( assembleWeakFormInDenseMode(assembler, context).release());
//	return shared_ptr<DiscreteBoundaryOperator<ResultType> >(DenseGlobalAssembler<BasisFunctionType, ResultType>:: assembleDetachedWeakForm(*this->dualToRange(), *this->domain(), *assembler, context).release());
		return shared_ptr<DiscreteBoundaryOperator<ResultType> >(DenseGlobalAssembler<BasisFunctionType, ResultType>:: assembleDetachedWeakFormPeter(str,*this->dualToRange(), *this->domain(), *assembler, context, solV, rhsV, wm).release());
//		return shared_ptr<DiscreteBoundaryOperator<ResultType> >(DenseGlobalAssembler<BasisFunctionType, ResultType>:: assembleDetachedWeakFormPeter(str,*this->dualToRange(), *this->domain(), *assembler, context).release());
//		return shared_ptr<DiscreteBoundaryOperator<ResultType> >(DenseGlobalAssembler<BasisFunctionType, ResultType>:: assembleDetachedWeakFormPeter(str,*this->dualToRange(), *this->domain(), assembler, context).release());
	}
	default:
		throw std::runtime_error("ACA or other not implemented");
	}
}

private:
    /** \cond PRIVATE */
    shared_ptr<CollectionOfKernels> m_kernels;
    shared_ptr<CollectionOfShapesetTransformations> m_testTransformations;
    shared_ptr<CollectionOfShapesetTransformations> m_trialTransformations;
    shared_ptr<TestKernelTrialIntegral> m_integral;
    /** \endcond */
};

} // namespace Bempp

#endif
