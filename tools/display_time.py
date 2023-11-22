import time
def display_time(t1, t2, message, display):
  ct = time.time()
  print(f"{message} {ct - t2:.2f} sec. \t Total: {ct - t1:.2f} sec.")
