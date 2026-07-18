module TMBGMagneticHybridWannier

include("TMBGMagneticHW.jl")
include("TMBGProjection.jl")
include("TMBGCommonBasis.jl")
include("TMBGMagneticGeometry.jl")
include("TMBGIntrinsicIdealGeometry.jl")
include("TMBGIdealComponent.jl")
include("TMBGSymmetricGaugeProjection.jl")
include("TMBGFigure4Projection.jl")

using .TMBGMagneticHW
using .TMBGProjection
using .TMBGCommonBasis
using .TMBGMagneticGeometry
using .TMBGIntrinsicIdealGeometry
using .TMBGIdealComponent
using .TMBGSymmetricGaugeProjection
using .TMBGFigure4Projection

export TMBGMagneticHW,
       TMBGProjection,
       TMBGCommonBasis,
       TMBGMagneticGeometry,
       TMBGIntrinsicIdealGeometry,
       TMBGIdealComponent,
       TMBGSymmetricGaugeProjection,
       TMBGFigure4Projection

end
