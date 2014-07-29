To use the PHP shared library, first check out the CoolProp sources using git:

```
git clone https://github.com/ibell/coolprop
```
move into the folder
```
cd coolprop/wrappers/PHP
```
Run the build script
```
./PHPBuilder
```
This copies the shared library into a location that PHP can load.  sudo is needed to make the copy

Modify the PHP.ini file that PHP will load to add
```
extension = "CoolProp.so" 
```
after the [PHP] header.

Now any php files should work

To enable some useful debugging at runtime (turn off for deployment), you can add

```
<?php
error_reporting  (E_ALL);
ini_set ('display_errors', true);
...
?>
```
to the top of the script.