package main.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class general_utils {
    // helper to stop the program with error and usage stuff
    public static void stop(String message) {
        System.err.println("Error: " + message);
        System.err.println("Usage: java -jar predictrjar [--probabilities] --model <model file> --format <txt|html> {--seq <fasta file>|--maf <multiple-alignment-folder>}");
        System.exit(1);
    }

    public static Path findScript() throws FileNotFoundException {
        Path current = Paths.get("").toAbsolutePath();

        Path propraRoot = null;
        Path walker = current;
        while (walker != null) {
            if (walker.getFileName() != null && walker.getFileName().toString().toLowerCase().contains("propra")) {
                propraRoot = walker;
                break;
            }
            walker = walker.getParent();
        }

        if (propraRoot == null) {
            throw new FileNotFoundException("propra directory could not be found!");
        }

        Path script = findFile(propraRoot, "discord_interaction", "create_ticket.py");
        if (script == null) {
            throw new FileNotFoundException("Could not find utils/create_ticket.py under: " + propraRoot);
        }

        return script;
    }

    public static Path findFile(Path root, String subdirectory, String filename) {
        File[] children = root.toFile().listFiles();
        if (children == null) return null;

        for (File child : children) {
            if (child.isDirectory()) {
                if (child.getName().equals(subdirectory)) {
                    File candidate = new File(child, filename);
                    if (candidate.exists()) return candidate.toPath();
                }

                Path result = findFile(child.toPath(), subdirectory, filename);
                if (result != null) return result;
            }
        }

        return null;
    }
}
