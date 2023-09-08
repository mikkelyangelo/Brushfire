import time

from field import Field

SIZE = 400
MAX_ITERS = 10000
START = (2, 2)
END = (98, 98)


def main():
    stat_time = time.time()
    field = Field(SIZE, MAX_ITERS, START, END)

    print("Parameters:")
    print("Size:", field.size, 'x', field.size)
    print("Step:", field.step)

    print("\nTime of creating object:", round((time.time() - stat_time), 2), 's')
    last_time = time.time()

    field.brush_fire(field.get_barriers())
    print("Time of BrushFire:", round((time.time() - last_time), 2), 's')
    last_time = time.time()

    field.gradient_descent()
    print("Time of gradient descent:", round((time.time() - last_time), 2), 's\n')
    last_time = time.time()

    field.best_first_search()
    print("\nTime of first search:", round((time.time() - last_time), 2), 's')

    print("Time of work:", round((time.time() - stat_time), 2), 's')
    # field.show_3d_capability()
    # field.show_2d_capability()


if __name__ == '__main__':
    main()
