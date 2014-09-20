import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

/**
 * @author Ethan Hart
 */
public class FIB {

    public static void main(String[] args) throws FileNotFoundException {
        Scanner scanner = new Scanner(new File(args[0]));
        String sNums = scanner.nextLine();

        String[] parts = sNums.split(" ");

        int months = Integer.parseInt(parts[0]);
        int rabs = Integer.parseInt(parts[1].trim());
        int r = rabbits(months, rabs);
        System.out.println(r);
    }

    private static int rabbits(int n, int k) {
        List<Integer> fib = new ArrayList<Integer>();
        for (int i=0;i<n;i++) {
            if (i < 2) {
                fib.add(1);
            }
            else {
                int adults = fib.get(fib.size() - 1); // Previous month's rabbits are adults
                int babies = fib.get(fib.size() - 2) * k; // Rabbits from 2 months ago have babies
                fib.add(adults + babies);
            }
        }

        return fib.get(fib.size() - 1);
    }
}
