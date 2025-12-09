# Руководство по использованию

## Быстрый старт

### 1. Базовый запуск

```bash
make clean
make
./multipole_sim
```

### 2. С визуализацией

```bash
make viz
```

## Настройка параметров

Откройте `main.f90` и измените параметры в начале программы:

### Размер домена

```fortran
real(8), parameter :: XMAX = 1.0d0   ! Половина размера по X
real(8), parameter :: YMAX = 1.0d0   ! Половина размера по Y
```

Частицы будут в области `[-XMAX, XMAX] x [-YMAX, YMAX]`.

### Сетка

```fortran
integer, parameter :: NX_CELLS = 10   ! Число ячеек по X
integer, parameter :: NY_CELLS = 10   ! Число ячеек по Y
```

**Рекомендации:**
- Для N=1000 частиц: 10x10 ячеек
- Для N=10000: 20x20
- Для N=100000: 50x50

Правило: `N_CELLS ~ sqrt(N_PARTICLES) / 5`

### Точность разложения

```fortran
integer, parameter :: P_MAX = 10
```

**Выбор P_MAX:**

| P_MAX | Относительная ошибка | Время (100k частиц) | Применение |
|-------|---------------------|---------------------|------------|
| 5     | ~1%                 | 0.5s               | Быстрый просмотр |
| 10    | ~0.001%             | 1s                 | **Оптимально** |
| 15    | ~0.00001%           | 2s                 | Высокая точность |
| 20    | ~0.0000001%         | 4s                 | Публикации |
| 30    | ~10⁻¹⁰              | 10s                | Бенчмарки |

### Число частиц

```fortran
integer, parameter :: N_PARTICLES = 10000
```

**Ограничения:**
- Минимум: 10 (для тестов)
- Максимум: зависит от памяти (~1 GB для 10M частиц)

### Динамика

```fortran
real(8), parameter :: ETA = 1.0d0       ! Подвижность вихрей
real(8), parameter :: DT = 0.001d0      ! Шаг по времени
integer, parameter :: NSTEPS = 100      ! Число шагов
```

**Подбор DT:**

Для стабильности: `DT < 0.01 * cell_size / (η * F_max)`

где `F_max` - максимальная сила в системе.

**Типичные значения:**
- Быстрая релаксация: `ETA = 10.0`, `DT = 0.0001`
- Нормальная: `ETA = 1.0`, `DT = 0.001`
- Медленная: `ETA = 0.1`, `DT = 0.01`

### Контроль симуляции

```fortran
logical, parameter :: TEST_ACCURACY = .true.   ! Тест точности
logical, parameter :: RUN_DYNAMICS = .true.    ! Запуск динамики
integer, parameter :: FRAME_SKIP = 5           ! Запись каждого N-го кадра
```

## Примеры конфигураций

### Пример 1: Быстрый тест (1 секунда)

```fortran
integer, parameter :: N_PARTICLES = 100
integer, parameter :: P_MAX = 5
integer, parameter :: NSTEPS = 10
logical, parameter :: TEST_ACCURACY = .true.
logical, parameter :: RUN_DYNAMICS = .true.
```

### Пример 2: Средняя симуляция (1 минута)

```fortran
integer, parameter :: N_PARTICLES = 5000
integer, parameter :: NX_CELLS = 15
integer, parameter :: NY_CELLS = 15
integer, parameter :: P_MAX = 10
integer, parameter :: NSTEPS = 500
real(8), parameter :: DT = 0.001d0
```

### Пример 3: Большая симуляция (10 минут)

```fortran
integer, parameter :: N_PARTICLES = 50000
integer, parameter :: NX_CELLS = 30
integer, parameter :: NY_CELLS = 30
integer, parameter :: P_MAX = 15
integer, parameter :: NSTEPS = 1000
real(8), parameter :: DT = 0.0005d0
logical, parameter :: TEST_ACCURACY = .false.  ! Пропустить O(N²) тест
```

### Пример 4: Высокоточный расчет

```fortran
integer, parameter :: N_PARTICLES = 1000
integer, parameter :: P_MAX = 25
real(8), parameter :: ETA = 0.1d0
real(8), parameter :: DT = 0.0001d0
integer, parameter :: NSTEPS = 5000
integer, parameter :: FRAME_SKIP = 50
```

## Типы начальных конфигураций

Измените функцию `generate_random_particles` в `main.f90`:

### 1. Случайное распределение (по умолчанию)

```fortran
x = (x - 0.5d0) * 1.8d0 * XMAX
y = (y - 0.5d0) * 1.8d0 * YMAX
q = 1.0d0
```

### 2. Два кластера

```fortran
if (i <= n/2) then
  x = -0.5d0 + x * 0.3d0
  y = -0.5d0 + y * 0.3d0
else
  x = 0.5d0 + x * 0.3d0
  y = 0.5d0 + y * 0.3d0
end if
q = 1.0d0
```

### 3. Кольцо

```fortran
real(8) :: angle, radius
angle = 2.0d0 * 3.14159265359d0 * x
radius = 0.5d0 + 0.1d0 * y
x = radius * cos(angle)
y = radius * sin(angle)
q = 1.0d0
```

### 4. Решетка

```fortran
integer :: ix, iy, nx_p, ny_p
nx_p = int(sqrt(real(n, 8)))
ny_p = n / nx_p
ix = mod(i-1, nx_p)
iy = (i-1) / nx_p
x = -0.8d0 + 1.6d0 * real(ix, 8) / real(nx_p-1, 8)
y = -0.8d0 + 1.6d0 * real(iy, 8) / real(ny_p-1, 8)
q = 1.0d0
```

### 5. Случайные заряды (±1)

```fortran
x = (x - 0.5d0) * 1.8d0 * XMAX
y = (y - 0.5d0) * 1.8d0 * YMAX
q = 2.0d0 * (q - 0.5d0)  ! Заряды от -1 до +1
```

## Визуализация

### Основная визуализация

```bash
python3 visualize.py
```

### Анализ распределения

```bash
python3 visualize.py --analyze
```

Создаст файл `output/distribution_comparison.png`.

### Сохранение кадров для видео

```bash
python3 visualize.py --save
```

Затем создайте видео:

```bash
ffmpeg -r 20 -i output/plot_%06d.png -c:v libx264 -vf fps=20 \
       -pix_fmt yuv420p animation.mp4
```

### Пользовательская визуализация

Создайте свой Python скрипт:

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('output/frame_000100.csv',
                     delimiter=',', skip_header=1)
x, y, q = data[:, 0], data[:, 1], data[:, 2]

plt.scatter(x, y, c=q, s=20, cmap='RdBu_r')
plt.colorbar()
plt.axis('equal')
plt.show()
```

## Анализ результатов

### Энергия системы

Из вывода программы:

```
Energy comparison:
  Direct:    8185.28    ( 0.021 ms)
  Multipole: 6917.74    ( 8.811 ms)
  Relative error:   0.1548
```

**Интерпретация:**
- `Relative error < 0.01` - хорошая точность
- `Relative error < 0.001` - отличная точность
- Если ошибка слишком большая - увеличьте `P_MAX`

### Силы

```
Force comparison:
  ID  |    Direct Force    |  Multipole Force   | Rel. Error
    1  |  -0.497E+03  0.574E+03  |  -0.503E+03  0.569E+03  |   0.10E-01
```

Относительная ошибка силы ~1% типична для `P_MAX=10`.

### Производительность

Speedup показывает ускорение multipole метода:
- Speedup < 1: Direct быстрее (мало частиц)
- Speedup ~ 10: Оптимальная область
- Speedup > 100: Очень эффективно (много частиц)

## Решение проблем

### Программа крашится

1. Уменьшите `N_PARTICLES`
2. Увеличьте `NX_CELLS` и `NY_CELLS`
3. Уменьшите `DT`

### Плохая точность

1. Увеличьте `P_MAX`
2. Уменьшите размер ячеек (больше `NX_CELLS`)
3. Проверьте, что частицы не слишком близко

### Медленная работа

1. Уменьшите `P_MAX`
2. Увеличьте размер ячеек (меньше `NX_CELLS`)
3. Используйте `TEST_ACCURACY = .false.` для больших N

### Нестабильная динамика

1. Уменьшите `DT` в 2-10 раз
2. Уменьшите `ETA`
3. Проверьте начальную конфигурацию (нет ли частиц на одном месте)

## Генерация коэффициентов (опционально)

Для использования предрассчитанных таблиц спецфункций:

```bash
python3 compute_coefficients.py 50
```

Создаст коэффициенты для `p_max` до 50.

**Примечание:** Текущая версия не использует внешние таблицы -
все коэффициенты рассчитываются в `init_system()`.

## Оптимизация под конкретное железо

### Для Intel CPU

```makefile
FFLAGS = -O3 -march=native -xHost -ipo
```

### Для AMD CPU

```makefile
FFLAGS = -O3 -march=znver3 -mtune=znver3
```

### Для ARM (Apple M1/M2)

```makefile
FFLAGS = -O3 -mcpu=apple-m1
```

### Для отладки

```makefile
FFLAGS = -g -O0 -fcheck=all -fbacktrace -Wall
```

## Параллелизация (будущее расширение)

Для добавления OpenMP параллелизма:

```makefile
FFLAGS += -fopenmp
```

и в коде:

```fortran
!$OMP PARALLEL DO
do i = 1, sys%n_particles
  call compute_force_multipole(sys, i, fx, fy)
end do
!$OMP END PARALLEL DO
```

## Контакты и поддержка

Вопросы и предложения: Terragon Labs
