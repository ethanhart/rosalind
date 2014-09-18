import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * @author Ethan Hart
 */
public class RNA {

    public static void main(String[] args) throws FileNotFoundException {
        String dna;
        Scanner scanner = new Scanner(new File(args[0]));
        dna = scanner.next();

        String rna = dna.replace("T", "U");

        System.out.printf(rna + "\n");
    }

}