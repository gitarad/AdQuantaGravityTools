import sys
import time

def log_with_timestamp(output_file):
    # Open the output file in append mode
    with open(output_file, 'a') as f:
        while True:
            # Read a line of output from stdin (program's output)
            line = sys.stdin.readline()
            
            # Get the current timestamp
            timestamp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            
            # If there's no more input, stop
            if not line:
                f.write(f"finished - {timestamp}")
                break
            
            
            # Prefix the line with the timestamp
            timestamped_line = f"{timestamp} - {line}"
            
            # Print the line to stdout
            sys.stdout.write(line)
            
            # Write the timestamped line to the file
            f.write(timestamped_line)
            f.flush()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python log_with_timestamp.py <output_file>")
        sys.exit(1)
    
    # The output file to save the logs to
    output_file = sys.argv[1]
    
    log_with_timestamp(output_file)