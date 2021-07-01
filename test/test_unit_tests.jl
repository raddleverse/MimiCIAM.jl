using MimiCIAM
using Test

@testset begin "Unit Testing"

    # test primary functionality
    m = MimiCIAM.get_model()
    run(m)

    # TODO test more functionalities, keyword args, etc.

end