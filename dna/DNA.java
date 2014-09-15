import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * @author Ethan Hart
 */
public class DNA {

    public static void main(String[] args) throws FileNotFoundException {
        String content;
        Scanner scanner = new Scanner(new File(args[0]));
        content = scanner.next();

        int a_count = content.length() - content.replace("A", "").length();
        int c_count = content.length() - content.replace("C", "").length();
        int g_count = content.length() - content.replace("G", "").length();
        int t_count = content.length() - content.replace("T", "").length();

        System.out.printf("%d %d %d %d\n", a_count, c_count, g_count, t_count);
    }

}
