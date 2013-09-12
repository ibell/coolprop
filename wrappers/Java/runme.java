// runme.java

public class runme {
    static {
        System.loadLibrary("CoolProp");
    }
    
    public static void main(String argv[]){
        System.out.println(CoolProp.Props("P",'T',300,'Q',0,"R134a"));
    }
}