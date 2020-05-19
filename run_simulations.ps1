

$dir =  "D:\Drive\Projetos\Doutorado2C\ValveEstimation\ValveEstimation\x64\Release\"
$out_dir = "D:\Drive\Projetos\Doutorado2C\test_data\"

$excitation = "sinusoidal"
$models = "karnopp-lugre-gms"

$params = " --directory=" + $out_dir + " --type=simulation --models=" + $models
$params += " --excitation=" + $excitation

&($dir + "ValveEstimation.exe") $params.Split(" ")