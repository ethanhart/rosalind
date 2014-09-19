import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
/**
 * @author Ethan Hart
 */
public class REVC {

    public static void main(String[] args) throws FileNotFoundException {
        String dna;
        String c = "";
        Scanner scanner = new Scanner(new File(args[0]));
        dna = scanner.next();
        //dna = "AAAACCCGGT";
        for (int i=0; i<dna.length(); i++) {
            char ch = dna.charAt(i);

            if (Character.toString(ch).equals("A")) {
                c += "T";
            }
            else if (Character.toString(ch).equals("T")) {
                c += "A";
            }
            else if (Character.toString(ch).equals("G")) {
                c += "C";
            }
            else if (Character.toString(ch).equals("C")) {
                c += "G";
            }
        }

        String revc = new StringBuffer(c).reverse().toString();
        System.out.println(revc);

    }

}
