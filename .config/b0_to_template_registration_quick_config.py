# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# AntsRegistration(MagicMonkeyBaseApplication) configuration
# -----------------------------------------------------------------------------

# Application traits configuration

c.AntsRegistration.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.AntsRegistration.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.AntsRegistration.log_level = 30

c.AntsRegistration.base_config_file = ""


# -----------------------------------------------------------------------------
# AntsConfiguration(MagicMonkeyConfigurable) configuration
# -----------------------------------------------------------------------------

c.AntsConfiguration.accross_modalities = True

c.AntsConfiguration.dimension = 3

c.AntsConfiguration.init_transform = [0, 1, 1]

c.AntsConfiguration.inlier_range = [0.005, 0.995]

c.AntsConfiguration.interpolation = "Linear"

c.AntsConfiguration.klass = "magic_monkey.config.ants.AntsConfiguration"

c.AntsConfiguration.match_histogram = False

c.AntsConfiguration.passes = [{
    "conv_eps": 1e-5,
    "conv_max_iter": [400, 200, 100],
    "conv_win": 30,
    "grad_step": 0.2,
    "klass": "magic_monkey.traits.ants.AntsRigid",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.2,
                64,
                "Regular",
                0.5
            ],
            "klass": "magic_monkey.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.8,
                64,
                "Regular",
                0.7
            ],
            "klass": "magic_monkey.traits.ants.MetricMI"
        }
    ],
    "shrinks": [
        12,
        8,
        4
    ],
    "smoothing": [
        6,
        4,
        2
    ]
}, {
    "conv_eps": 1e-6,
    "conv_max_iter": [500, 300, 200, 100],
    "conv_win": 20,
    "grad_step": 0.1,
    "klass": "magic_monkey.traits.ants.AntsAffine",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.2,
                64,
                "Regular",
                0.5
            ],
            "klass": "magic_monkey.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.8,
                64,
                "Regular",
                0.7
            ],
            "klass": "magic_monkey.traits.ants.MetricMI"
        }
    ],
    "shrinks": [
        12,
        8,
        4,
        2
    ],
    "smoothing": [
        6,
        4,
        2,
        1
    ]
}, {
    "conv_eps": 1e-6,
    "conv_max_iter": [200, 200, 100, 100, 30],
    "conv_win": 10,
    "grad_step": 0.1,
    "var_penality": 3,
    "var_total": 0,
    "klass": "magic_monkey.traits.ants.AntsSyN",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.1,
                64,
                "Regular",
                0.5
            ],
            "klass": "magic_monkey.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.6,
                64,
                "Regular",
                0.7
            ],
            "klass": "magic_monkey.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 2,
            "args": [
                0.3,
                64,
                "Regular",
                0.7
            ],
            "klass": "magic_monkey.traits.ants.MetricMI"
        }
    ],
    "shrinks": [
        12,
        8,
        4,
        2,
        1
    ],
    "smoothing": [
        6,
        4,
        2,
        1,
        0
    ]
}, {
    "conv_eps": 1e-6,
    "conv_max_iter": [40, 20],
    "conv_win": 10,
    "grad_step": 0.1,
    "var_penality": 3,
    "var_total": 0,
    "klass": "magic_monkey.traits.ants.AntsSyN",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.1,
                8,
                "Regular",
                0.5
            ],
            "klass": "magic_monkey.traits.ants.MetricCC"
        },
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.6,
                8,
                "Regular",
                0.7
            ],
            "klass": "magic_monkey.traits.ants.MetricCC"
        },
        {
            "target_index": 0,
            "moving_index": 2,
            "args": [
                0.3,
                8,
                "Regular",
                0.7
            ],
            "klass": "magic_monkey.traits.ants.MetricCC"
        }
    ],
    "shrinks": [
        2,
        1
    ],
    "smoothing": [
        1,
        0
    ]
}]

c.AntsConfiguration.use_float = False
