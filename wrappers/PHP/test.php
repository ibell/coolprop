<?php
require "CoolProp.php";

$p = 101325;
$Q = 1.0;
$T = PropsSI("T","P",$p,"Q",$Q,"Water");
print "NBP of water is $T\n";
?>
