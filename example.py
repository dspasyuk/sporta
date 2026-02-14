import sporta
import pandas as pd

# Example: Process a single file
try:
    # Process the first frame of a master file
    result = sporta.process_file("./data/insul_1_master.h5", index=1, threshold=6)
    print("Analysis Result:")
    print(result)
except Exception as e:
    print(f"Error: {e}")

# Example: Process multiple frames
results = []
for i in range(1, 4):
    try:
        res = sporta.process_file("./data/insul_1_master.h5", index=i)
        res['Frame'] = i
        results.append(res)
    except Exception as e:
        print(f"Skipping frame {i}: {e}")

if results:
    df = pd.DataFrame(results)
    print("\nBatch Results:")
    print(df)
