import json
from math import exp

import numpy as np
import matplotlib.pyplot as plt
from treelib import Tree
import sys

from cell import Cell

np.set_printoptions(threshold=sys.maxsize)
PRIMES = (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 43, 61, 67, 71, 73)


class Field:
    """Класс поля с препятствиями и нахождения пути для прохода по нему."""

    def __init__(self, size=100, iters=2000, start=(1, 1), end=(99, 99)):
        """Инициализация объекта поля."""

        self.size = size

        self.start = self._prepare_start_end(start)
        self.end = self._prepare_start_end(end)

        self.iters = iters
        self.step = 100.0 / size
        self.field = np.zeros((size, size), dtype=Cell)
        self._fill_field()
        self.beta = 1

    def _prepare_start_end(self, obj):
        """Преобразование координат """

        return int(obj[0] * (self.size / 100) - 0.01), int(abs(self.size - obj[1] * (self.size / 100)) - 0.01)

    def _fill_field(self):
        """Заполнение массива клеток поля объектами класса Клетки."""

        for y in range(self.size):
            for x in range(self.size):
                self.field[y][x] = Cell(x, y, self.end)

    def _stepping(self, x):
        """Получение числа кратного шагу."""

        return round((x - x % self.step), 2)

    def get_barriers(self):
        """Расстановка препятствий на поле."""

        with open('ATP/sample_output_file1.json') as f:
            barriers = json.load(f)

        return_list = []
        # Словарь списков областей каждого препятствия
        list_of_ones = {str(i): [] for i in range(len(barriers)-3)}
        for i_pol, polygon in enumerate(barriers[3:]):
            points_list = []

            # Создание списка координат вершин фигуры
            for point in polygon['points']:
                x = point['x']
                y = point['y']

                if x == 100.0:
                    x = 100.0 - 0.001
                if y == 100.0:
                    y = 100.0 - 0.001

                points_list.append([x, y])
            points_list = np.array(points_list)

            # Расстановка вершин фигуры в поле
            for point in points_list:
                x = self._stepping(point[0])
                y = self._stepping(point[1])

                list_of_ones[str(i_pol)].append((y, x))

            # Прохождение по всем парам вершин фигуры (граням)
            # для обозначения вертикальных и наклонных граней
            for i in range(points_list.shape[0]):
                j = i + 1 - (i // (points_list.shape[0] - 1)) * points_list.shape[0]

                y1, y2 = points_list[i][1], points_list[j][1]
                x1, x2 = points_list[i][0], points_list[j][0]

                # Прохождение по каждой точке на прямой с шагом
                for y in np.arange(self._stepping(min((y1, y2))) + self.step, self._stepping(max((y1, y2))), self.step):
                    x = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                    x = self._stepping(x)

                    list_of_ones[str(i_pol)].append((y, x))

            # Заполнение горизонтальных линий фигур
            for y in np.arange(self._stepping(np.min(points_list[:, 1])),
                               self._stepping(np.max(points_list[:, 1])), self.step):
                counter = 0
                cur_line = []

                # Читаем кол-во границ на горизонтальной
                for yx in list_of_ones[str(i_pol)]:
                    if yx[0] == y:
                        counter += 1
                        cur_line.append(yx[1])

                # Если есть границы, заполняем пространство между ними
                if counter >= 2:
                    for x in np.arange(min(cur_line), max(cur_line), self.step):
                        if (y, x) not in list_of_ones[str(i_pol)]:
                            list_of_ones[str(i_pol)].append((y, x))

            # Наносим все по-отдельности фигуры на поле
            for point in list_of_ones.keys():
                for one in list_of_ones[point]:
                    x = int(one[1] / self.step)
                    y = int(one[0] / self.step)

                    self.field[self.size - y - 1][x].distance = 1
                    return_list.append((self.size - y - 1, x))

        return return_list

    def _get_neighbours(self, cell):
        """Получение соседей клетки."""

        x, y = cell.x, cell.y
        # neighs = ((x-1, y), (x, y+1), (x+1, y), (x, y-1))
        neighs = []
        for xs in range(x-1, x+2):
            for ys in range(y-1, y+2):
                if xs != x and ys != y:
                    neighs.append((xs, ys))

        return tuple(self.field[xy[1]][xy[0]] for xy in neighs
                     if xy[0] in range(self.size) and xy[1] in range(self.size))

    def brush_fire(self, list_of_ones):
        """Алгоритм Brushfire. Нахождение расстояния до ближайшей преграды дял каждой точки."""

        borders = []
        for xy in list_of_ones:
            borders.append(xy)

        # Расстановка расстояний для границ поля
        for y in range(self.size):
            for x in range(self.size):

                if x in (0, self.size-1) or y in (0, self.size-1):
                    self.field[y][x].distance = 1
                    borders.append((x, y))

        length = 1
        while True:
            length += 1

            new_borders = []

            for xy in borders:
                cell = self.field[xy[0]][xy[1]]

                for neigh in self._get_neighbours(cell):
                    if neigh.distance == 0:
                        neigh.distance = length
                        new_borders.append((neigh.y, neigh.x))

            if len(new_borders) == 0:
                break

            borders = new_borders

    def _get_beta(self):
        """Нахождение Beta для градиентного спуска."""

        max_length = np.max(self.field).distance
        for num in PRIMES:
            if max_length < num:
                break
            self.beta *= num

    def gradient_descent(self):
        """Градиентный спуск для нахождения 'потенциала' каждой точки."""

        self._get_beta()

        # Расчет сил притяжения для каждой точки
        for y in range(self.size):
            for x in range(self.size):
                cell = self.field[y][x]
                # cell.capability = int(self.beta / abs(0.5 - cell.distance))
                # cell.capability = 100000 * pow((1 / cell.distance), 2)
                cell.capability = 100 * pow((1 / (0.25 - cell.distance)), 2)
                # cell.capability = abs((1 / (0.99 - cell .distance)))
                # cell.capability = exp(25 * (1 / cell.distance))
                cell.capability += 1 * abs(self.end[0] - x) + abs(self.end[1] - y)

    @staticmethod
    def _find_min(arr):
        """Нахождение точки с минимальным 'потенциалом' из списка."""

        min_ = np.inf
        min_cell = np.inf
        # print(len(arr))
        for cell in arr:
            if cell.distance != 1:
                # print((cell.x, cell.y), cell.capability, min_)
                if cell.capability == min_ and cell.evcl_distance < min_cell.evcl_distance or cell.capability < min_:
                    min_ = cell.capability
                    min_cell = cell

        return min_cell

    def best_first_search(self):
        """Алгоритм поиска от наилучшего. Поиск пути до конечной точки."""

        start = self.field[self.start[1]][self.start[0]]

        # Клетки принадлежащие фигурам
        skipped = set()
        # Дерево пути
        all_tops = Tree()
        # Начальная точка
        all_tops.create_node(identifier=str((start.x, start.y)))

        neighs = self._get_neighbours(start)
        for neigh in neighs:
            try:
                all_tops.create_node(identifier=str((neigh.x, neigh.y)), parent=str((start.x, start.y)))
            except:
                pass

        quite = False
        for _ in range(self.iters):
            # while True:
            tops = set([self.field[eval(nod.tag)[1]][eval(nod.tag)[0]] for nod in all_tops.leaves()])

            top = self._find_min(tops - skipped)
            # Проверка на отсутствие пути
            if isinstance(top, float):
                print("No path!")
                break

            cnt = 0
            neighbours = self._get_neighbours(top)
            for neigh in neighbours:
                try:
                    all_tops.create_node(identifier=str((neigh.x, neigh.y)), parent=str((top.x, top.y)))
                except:
                    cnt += 1

                # Условие достижения финальной точки
                if neigh.x == self.end[0] and neigh.y == self.end[1]:
                    quite = True

            if quite:
                break

            # Добавление точки принадлежащей фигуре
            if cnt == len(neighbours):
                skipped.add(top)

        # all_tops.show()
        xy = self._find_min_dist(all_tops.all_nodes())
        print("End point:", xy)

        self._create_way(list(all_tops.rsearch(str(xy))))

    def _find_min_dist(self, leaves):
        """Нахождение точки с минимальным расстоянием до конечной точки."""

        dist = lambda x_, y_: abs(self.end[0] - x_) + abs(self.end[1] - y_)
        min_dist = self.size ** 2
        xy = (1, 198)

        for leaf in leaves:
            x = eval(leaf.tag)[0]
            y = eval(leaf.tag)[1]
            if dist(x, y) <= min_dist:
                min_dist = dist(x, y)
                xy = (x, y)

        return xy

    def _create_way(self, tops):
        """Преобразование в .json файл для передачи его в ATP генератор."""

        tops_json = [
            {
                "type": "polyline",
                "size": 1.0,
                "color": 4572354,
                "points": []
            }
        ]

        for point in tops:
            point = eval(point)
            tops_json[0]['points'].append({
                "type": "point",
                "x": float(point[0]) * (100 / self.size),
                "y": float(self.size - point[1] - 1) * (100 / self.size)
            })

        with open('ATP/sample_user_solve1.json', 'w') as f:
            json.dump(tops_json, f, indent=4, ensure_ascii=False)

    def show_3d_capability(self):
        # Creating dataset
        x = np.arange(0, self.size).reshape(1, -1)
        x = np.repeat(x, self.size, axis=0)
        y = np.arange(0, self.size).reshape(1, -1)
        y = np.repeat(y, self.size, axis=0).T
        z = np.zeros((self.size, self.size))

        for xi in range(self.size):
            for yi in range(self.size):
                z[xi][yi] = self.field[xi][yi].capability
        z = self._normalize(z)

        ax = plt.axes(projection='3d')
        ax.plot_surface(x, y, z)
        plt.show()

    def show_2d_capability(self):
        """Создание 2д карты."""

        z = np.zeros((self.size, self.size))

        for xi in range(self.size):
            for yi in range(self.size):
                z[xi][yi] = self.field[xi][yi].capability
                # if self.field[xi][yi].distance == 1:
                #     z[xi][yi] = 1
                # else:
                #     z[xi][yi] = 0

        # z = self._normalize(z)

        plt.imshow(z, cmap='inferno', vmin=np.min(z), vmax=np.max(z))   # np.max(z)
        plt.show()

    @staticmethod
    def _normalize(field):
        """Нормализация данных для отображения на графике."""

        field -= np.min(field)
        field -= np.mean(field)
        field /= np.max(field)

        return field
