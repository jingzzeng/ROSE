import copy, dcor, torch, scipy, os, time, _ctypes, warnings, matplotlib
import numpy as np
from numbers import Integral, Real
from torch import nn
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils import check_X_y, check_array, check_random_state
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.validation import check_is_fitted
from sklearn.preprocessing import OneHotEncoder
from matplotlib import cm




class NSDR(TransformerMixin, BaseEstimator):

    # convert NSDR to sklearn compatible form. GridCV, pipeline can be directly applied.
    # all check estimator test pass except datatype. float32 is assume here to accelerate computation.

    _parameter_constraints = {
    "initial_blambda": [Interval(Real, 0, None, closed="left"), StrOptions({'auto'})],
    "max_dim": [Interval(Integral, 0, None, closed="left")],
    "method": [StrOptions({"fnorm", "seq", "fnormsquare"})],
    "retrain": ["boolean"],
    "adaptive_cv": ["boolean"],
    "device": [StrOptions({'cpu', "cuda:0", "cuda:1"})], ## for later check should number to be specified
    "debug": ["boolean"],
    "iter_num": [Interval(Integral, 1, None, closed="left")],
    "early_stop": ["boolean"],
    "y_mode": [Interval(Integral, 0, 5, closed="both")],
    "batch_size": [Interval(Integral, 12, None, closed="left")],
    "lr": [Interval(Real, 0, 1, closed="both")],
    "blambda_tilde": [Interval(Real, 0, None, closed="left"), StrOptions({"auto"})],
    "random_state": ["random_state"],
    "optimizer": [StrOptions({"adam", "sgd"})],
    "categorical_y": ["boolean"],
    "y_mode":[Interval(Integral, 0, 5, closed="both")]
    }
    

    def __init__(self, neural_network=id(None), initial_blambda="auto", max_dim=1,
        method="seq", iter_num=100, early_stop=False, y_mode=0, random_state=None, 
        categorical_y=False, batch_size = 100, lr = 1e-3, blambda_tilde="auto", 
        adaptive_cv=False, retrain=True, optimizer="adam",device = 'cpu', 
        debug=False, debug_logfilename="default.txt"):
        self.set_params(neural_network = neural_network)
        self.initial_blambda = initial_blambda  # initial value of lambda parameter
        self.max_dim = max_dim
        self.retrain = retrain
        self.debug = debug
        self.adaptive_cv = adaptive_cv
        self.device = device
        self.iter_num = iter_num
        self.categorical_y = categorical_y
        self.early_stop = early_stop
        self.y_mode = y_mode
        self.batch_size = batch_size
        self.optimizer = optimizer
        self.random_state = random_state
        self.lr = lr
        self.blambda_tilde = blambda_tilde
        self.method = method
        if isinstance(debug_logfilename, str):
            self.debug_logfilename = debug_logfilename
        else:
            raise ValueError("debug_logfilename should be a string specifying the directory")

    def set_params(self, **params):
        for key, value in params.items():
            if key == "neural_network":
                if isinstance(value, nn.Module):
                    # store its id to pass estimator test of sklearn
                    self.neural_network = id(value)
                    self.neural_network_ = value
                elif isinstance(value, str):
                    try:
                        self.neural_network = id(eval(value))
                    except Exception as e:
                        raise ValueError(f"Evaluation fails for {value}")
                elif isinstance(value, int):
                    if value < int(hex(id(None))[:4] + "0" * (len(hex(id(None))) - 4), 16):
                        raise ValueError("u should not pass int as parameters")
                    # else it's address.
                    self.neural_network = value
                else:
                    raise ValueError("neural_network should be specified correctly or leave it as default")
            else:
                super().set_params(**{key: value})
        return self

    def __getstate__(self):
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        if self.debug:
            del state['file_']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
    

    def _output_to_debugfile(self, word):
        if self.debug:
            print(word)
            self.file_.write(word)
            self.file_.write("\n")
            self.file_.flush()


    def _gmdd_dist(self, y, y_mode = 0):
        temp_y = y.reshape(y.shape[0], -1)
        tmp_dist_y = torch.cdist(temp_y, temp_y)
        if y_mode == 0:
            return tmp_dist_y
        elif y_mode == 1:
            return tmp_dist_y ** 2 / 2 # 6 th case in paper
        elif y_mode == 2:
            return tmp_dist_y ** 2 / (1 + tmp_dist_y ** 2) # 3rd case in paper with lambda = 1
        elif y_mode == 3:
            return 1 - torch.exp(- tmp_dist_y ** 2 / 2) # 4th case
        elif y_mode == 4:
            return torch.log(1 + tmp_dist_y ** 2 / 2) # 9th case
        elif y_mode == 5:
            return tmp_dist_y ** 1.5 # 7 th case with alpha = 1.5
        elif y_mode == 6:
            ## distance for categorical y
            return (tmp_dist_y != 0).to(torch.float32)


    def load_pretrained_model(self, path_list):
        self.preloadnet_list_ = []
        if self.neural_network == id(None):
            raise Exception("u should specify the neural_network paramter before load the pretrained model")
        for path in path_list:
            net = copy.deepcopy(_ctypes.PyObj_FromPtr(self.neural_network))
            net = net.to(self.device)
            net.load_state_dict(torch.load(path))
            self.preloadnet_list_.append(net)


    def score(self, X, y):
        ypred = self.transform(X)
        return dcor.distance_correlation(ypred, y)


    def _gmdd_value(self, yhat, y):
        ## return the gmdd value based on the unbiased version of E[(yhat - E yhat)(\tilde{yhat} - E yhat) \|y-\tilde{y}\|] 
        n_samples = y.shape[0]
        H = torch.eye(n_samples) - torch.ones(n_samples, 1) @ torch.ones(1, n_samples) / n_samples
        H = H.to(self.device)
        tmp_dist_y = self._gmdd_dist(y, y_mode=self.y_mode)
        tmp_yhat_mean = torch.mean(yhat, axis = 0).reshape(1, -1)
        tmp_centered_yhat = yhat - tmp_yhat_mean.repeat(n_samples, 1)
        A = H @ tmp_centered_yhat @ tmp_centered_yhat.T @ H
        B = H @ tmp_dist_y @ H
        gmdd = torch.mean(A * B)
        return gmdd


    @staticmethod
    def generate_default_net(input_dim, first_hidden_width=None, last_hidden_width=16, increase_count=2,  output_dim=1):
        if first_hidden_width is None:
            first_hidden_width = 2 ** (int(np.log2(input_dim)) + 1)
        net = nn.Sequential()
        net.append(nn.Linear(input_dim, first_hidden_width))
        net.append(nn.ReLU())
        next_width = first_hidden_width
        for i in range(increase_count - 1):
            net.append(nn.Linear(next_width, next_width * 2))
            net.append(nn.ReLU())
            next_width = next_width * 2
        for i in range(increase_count - 1):
            net.append(nn.Linear(next_width, next_width // 2))
            net.append(nn.ReLU())
            next_width = next_width // 2
        if last_hidden_width is not None:
            while(next_width // 2 >= last_hidden_width):
                net.append(nn.Linear(next_width, next_width // 2))
                net.append(nn.ReLU())
                next_width = next_width // 2
        net.append(nn.Linear(next_width, output_dim))
        return net

    @staticmethod
    def init_xv_uniform(m):
        if type(m) == nn.Linear:
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)

    def _generate_pred_inbatch(self, net, data_loader):
        net.eval()
        pred = torch.Tensor([]).to(self.device)
        with torch.no_grad():
            for data in data_loader:
                # data loader can yield (x, y) tuple or x tuple 
                temp_x = data[0]
                yhat = net(temp_x)
                pred = torch.cat((pred, yhat), 0)
        return pred



    def _dist_cor(self, neural_network, data_loader):
        # return average of distance correlation of each batch
        neural_network.eval()
        dc = 0
        with torch.no_grad():
            for temp_x, temp_y in data_loader:
                yhat = neural_network(temp_x)
                dc += dcor.distance_correlation(yhat.cpu().numpy(), temp_y.cpu().numpy())
        dc = dc / len(data_loader)
        return dc



    def _gmdd_of_loader(self, neural_network, data_loader, standardize=True, batch_average=False):
        neural_network.eval()
        loader_gmdd = 0
        with torch.no_grad():
            if data_loader.batch_size < 10:
                raise Exception("batch size of data_loader in _gmdd_of_loader should be at least 10")
            if not batch_average and len(data_loader.dataset) > 5000:
                raise Exception("sample size of _gmdd_of_loader should be at most 5k if calculated without batch_average")
            if batch_average:
                with torch.no_grad():
                    for temp_x, temp_y in data_loader:
                        yhat = neural_network(temp_x)
                        if standardize:
                            loader_gmdd += self._gmdd_value(yhat, temp_y) / torch.sum(torch.diag(torch.cov(yhat.T).reshape(yhat.shape[1], -1))) # reshape to avoid bug when yhat is N * 1 matrix
                        else:
                            loader_gmdd += self._gmdd_value(yhat, temp_y)
                loader_gmdd = loader_gmdd / len(data_loader)
            else:
                yhat = []
                y = []
                for temp_x, temp_y in data_loader:
                    # Assuming each batch contains a tuple (inputs, labels)
                    y.append(temp_y)
                    yhat.append(neural_network(temp_x))

                # Concatenate the data into a single tensor or list
                y = torch.cat(y, dim=0)
                yhat = torch.cat(yhat, dim=0)
                if standardize:
                    loader_gmdd += self._gmdd_value(yhat, y) / torch.sum(torch.diag(torch.cov(yhat.T).reshape(yhat.shape[1], -1))) # reshape to avoid bug when yhat is N * 1 matrix
                else:
                    loader_gmdd += self._gmdd_value(yhat, y)

        return loader_gmdd

    def get_sample_gmdd_value(self, X, y, standardize=True, batch_average=None):
        check_is_fitted(self)
        if batch_average is None:
            if X.shape[0] <= 5000:
                batch_average = False
            else:
                batch_average = True
        X, y = self._validate_data(X, y, allow_nd=True,  multi_output=True,  y_numeric=True, reset=False)
        X = torch.Tensor(X).to(self.device)
        y = torch.Tensor(y).to(self.device) # conversion from numpy to torch
        if self.categorical_y and (y.ndim == 1 or y.shape[1] == 1):
            y = torch.nn.functional.one_hot(y.to(torch.long))
            y = y.to(torch.float32).reshape(y.shape[0], -1)
        gmdd = 0
        data_loader = torch.utils.data.DataLoader(torch.utils.data.TensorDataset(X, y), batch_size=self.batch_size, shuffle=True)  ## not reshuffle after each load
        if self.method == "seq":
            for i in range(self.max_dim):
                gmdd += self._gmdd_of_loader(self.estimators_[i], data_loader, standardize, batch_average)
        else:
            gmdd = self._gmdd_of_loader(self.estimators_[0], data_loader, standardize, batch_average)
        return gmdd.data.cpu().numpy()

    def get_sample_gmdd_matrix(self, X, y):
        ## return 
        check_is_fitted(self)
        X, y = self._validate_data(X, y, allow_nd=True, multi_output=True, y_numeric=True, reset=False)
        y = torch.Tensor(y).to(self.device)
        if self.categorical_y and (y.ndim == 1 or y.shape[1] == 1):
            y = torch.nn.functional.one_hot(y.to(torch.long))
            y = y.to(torch.float32).reshape(y.shape[0], -1)

        n_samples = y.shape[0]
        yhat = torch.Tensor(self.transform(X)).to(self.device)
        if n_samples > 5000:
            warnings.warn("It's not preferred to calculate the empirical gmdd matrix with large sample size")
        tmp_yhat_mean = torch.mean(yhat, axis = 0).reshape(1, -1)
        tmp_centered_yhat = yhat - tmp_yhat_mean.repeat(n_samples, 1)
        gmddsize = yhat.shape[1]
        gmdd = torch.zeros((gmddsize, gmddsize))
        cdisty = self._gmdd_dist(y, y_mode=self.y_mode)
        for i in range(gmddsize):
            for j in range(i, gmddsize):
                gmdd[i, j] = torch.mean(torch.outer(tmp_centered_yhat[:, i], tmp_centered_yhat[:, j]) * cdisty)
        tindx = torch.tril_indices(gmddsize, gmddsize, -1)
        gmdd[tindx[0], tindx[1]] = gmdd.T[tindx[0], tindx[1]] ## convert upper tri to sysm
        return gmdd.data.cpu().numpy()

    @staticmethod
    def torch_normsq(x):
        return (torch.norm(x)) ** 2

    def fit(self, X, y):
        # validity checker
        torch.set_default_dtype(torch.float32)
        self._validate_params() # sanity check of initial params
        self._rng = check_random_state(self.random_state) # get numpy random state. use it to generate consistent random number
        torch.manual_seed(self._rng.randint(0, 0xffff_ffff))  # set global seed of torch. self._rng.randint() return the same value if random state set equal
        # X, y = check_X_y(X, y, allow_nd=True, multi_output=True, y_numeric=True)
        X, y = self._validate_data(
            X, y, allow_nd=True, y_numeric=True, multi_output=True, ensure_min_samples=10
        )

        if self.categorical_y and np.unique(y, axis=0).shape[0] > y.shape[0] / 2: # each class get average two sample
            raise ValueError(f"y has {y.shape[0]} samples and {np.unique(y, axis=0)} category. y should be numeric rather than categorical")
        elif not self.categorical_y and y.shape[0] / np.unique(y, axis=0).shape[0] > 100:
            raise ValueError(f"y may be categorical. u specify categorical_y as True to get better performance.")

        if self.categorical_y:
            self.y_mode = 6

        if self.blambda_tilde == "auto":
            # set initial blambda_tilde_ and reset it after first direction fitted. 
            self.blambda_tilde_ = 1 # the actual blambda tilde
        else:
            self.blambda_tilde_  = self.blambda_tilde

        if self.neural_network == id(None):
            if self.method == "seq":
                self.nn_output_dim_ = 1
            else:
                self.nn_output_dim_ = y.shape[1]
            self.neural_network_ = self.generate_default_net(X.shape[1], 2 ** (int(np.log2(X.shape[1]) +1)), 16, 2, self.nn_output_dim_)
        else:
            self.neural_network_ = _ctypes.PyObj_FromPtr(self.neural_network) ## get the neural network from its id
        keyrev = reversed(self.neural_network_.state_dict())
        self.nn_output_dim_ = self.neural_network_.state_dict()[next(keyrev)].shape[0]

        if self.initial_blambda == "auto":
            if self.nn_output_dim_ == 1:
                # set initial_blambda_ small if we fit nn with output dim 1
                self.initial_blambda_ = 0.1 ** 6
            else:
                self.initial_blambda_ = 1 ## initial blambda of fnorm should be the estimate of   $\sqrt{d} lambda_1^*$. rather than auto
        else:
            self.initial_blambda_ = self.initial_blambda

        if self.method == "seq" and self.nn_output_dim_ > 1:
            raise ValueError("Output dimension of neural network is greater than 1, u should use fnorm method instead of sequential method")
        if self.max_dim > 1 and self.method != "seq":
            raise ValueError("Sequential method is preferred u train NSDR with max_dim greater than 1")
        if self.debug:
            if (not os.path.isfile(self.debug_logfilename)) and (not os.path.exists(os.path.dirname(self.debug_logfilename))):
                raise Exception("debug logfile not exist")
            self.file_ = open(self.debug_logfilename, "w")
            self.file_.write(f"NSDR start at {time.ctime(time.time())} \n")
            self.file_.flush()

        if self.method == "fnorm":
            self.normfunc_ = torch.norm
            self.lower_cov_ = 0
            self.upper_cov_ = np.inf
        elif self.method == "seq":
            self.normfunc_ = self.torch_normsq
            self.lower_cov_ = 0
            self.upper_cov_ = np.inf  #maximum cov empirical estimate can have
        elif self.method == "fnormsquare":
            self.normfunc_ = self.torch_normsq
            self.lower_cov_ = 0
            self.upper_cov_ = np.inf
        else:
            raise Exception("norm should be fnorm or seq")

        # actual fit method

        if self.retrain:
            if not hasattr(self, 'preloadnet_list_'):
                self.preloadnet_list_ = []
        else:
            if not hasattr(self, 'preloadnet_list_'):
                raise Exception("u should preload paramter of net list before fitting")
        self.estimators_ = []
        X = torch.Tensor(X).to(self.device)
        y = torch.Tensor(y).to(self.device) # conversion from numpy to torch
        if self.categorical_y and (y.ndim == 1 or y.shape[1] == 1):
            y = torch.nn.functional.one_hot(y.to(torch.long))
            y = y.to(torch.float32).reshape(y.shape[0], -1)


        self._fit_stages(X, y)

        return self




    def _fit_stages(self, X, y):
        self.estimators_.extend(self.preloadnet_list_)
        shuffled_dataset = torch.utils.data.Subset(torch.utils.data.TensorDataset(X, y), torch.randperm(len(X)).tolist())
        train_loader = torch.utils.data.DataLoader(shuffled_dataset, batch_size=self.batch_size, shuffle=False)
        for i in range(len(self.preloadnet_list_), self.max_dim):
            self.estimators_.append(copy.deepcopy(self.neural_network_)) # set net to list

        for i in range(self.max_dim):
            self.estimators_[i] = self.estimators_[i].to(self.device)

        for i in range(self.max_dim):
            self._output_to_debugfile(f"{i}th net is fitting")
            if i >= len(self.preloadnet_list_) or (i < len(self.preloadnet_list_) and self.retrain):
                # refit if it's a initialized net or specify it.
                fitted = self._fit_stage(i, X, y)
            else:
                fitted = self.estimators_[i]
            fitted.eval() # set it to evaluation mode.
            pred = self._generate_pred_inbatch(fitted, train_loader)
            if self.method == "seq":
            # if True:
                # normalize for seq method. 
                # Fnorm for net with output dimension may have collinearity and thus degeneration of covariance 
                # And inverse of matrix with zero eigvalue will incur large deviation. 
                keyrev = reversed(fitted.state_dict())
                fitted_lastlayer_bias = fitted.state_dict()[next(keyrev)]
                fitted_lastlayer_weight = fitted.state_dict()[next(keyrev)]
                fitted_lastlayer_weight.copy_(fitted_lastlayer_weight/torch.sqrt(torch.cov(pred.T)))
                pred = self._generate_pred_inbatch(fitted, train_loader)  ## refit after weight have alter
                fitted_lastlayer_bias.copy_(fitted_lastlayer_bias - torch.mean(pred, axis=0))

            if self.method == "seq" and i == 0 and self.blambda_tilde == "auto":
                if len(train_loader.dataset) <= 5000:
                    batch_average = False
                else:
                    batch_average = True
                self.blambda_tilde_ = (-self._gmdd_of_loader(fitted, train_loader, standardize=True, batch_average=batch_average)).detach()
                self._output_to_debugfile(f"self.blambda_tilde_ {self.blambda_tilde_}")

            self.estimators_[i] = fitted
        pred = self._generate_pred_inbatch(fitted, train_loader)
        if pred.shape[1] > 1 and self.debug:
            self._output_to_debugfile(f"""corelation of net is {torch.corrcoef(pred.T)} and 
                the fnorm of corelation is {torch.norm(torch.corrcoef(pred.T) - torch.eye(pred.shape[1], device=self.device)) / pred.shape[1]} \n""")
                

    def _fit_stage(self, i, X, y):
        ## fit the i th net of estimators
        neural_network = self.estimators_[i]
        training_dataset = torch.utils.data.TensorDataset(X, y)
        train_loader = torch.utils.data.DataLoader(training_dataset, batch_size=self.batch_size, shuffle=True)
        tr_dataset, cv_dataset = torch.utils.data.random_split(training_dataset, [0.5, 0.5])
        tr_loader = torch.utils.data.DataLoader(tr_dataset, batch_size=self.batch_size, shuffle=True)
        cv_loader = torch.utils.data.DataLoader(cv_dataset, batch_size=self.batch_size, shuffle=True)  # used for evaluate distance correlation, large batch for stability. 
        precision = 6
        upper_blambda = 10 ** precision
        current_blambda = self.initial_blambda_
        lower_blambda = 0.1 ** precision
        iterated = set([current_blambda])
        old_discor = 0
        dec_count = 0
        best_blambda = current_blambda
        best_dis = 0
        while self.adaptive_cv:
            ## tuning process of blambda
            self.blambda_ = current_blambda
            fitted_nn = self._fit_one_neural_network(i, tr_loader)
            new_discor = self._dist_cor(fitted_nn, cv_loader)
            self._output_to_debugfile(f"blambda is {current_blambda} discor is {new_discor} cov sum is {self.cov_sum_[-1]}")
            if new_discor > best_dis and self.lower_cov_ <= self.cov_sum_[-1] <= self.upper_cov_:
                best_dis = new_discor
                best_blambda = current_blambda
                self._output_to_debugfile(f"best discor is {best_dis} best blambda is {best_blambda}")
            if self.cov_sum_[-1] < 10 ** 2: # if cov sum is large lambda should be increased
                upper_blambda = current_blambda
                current_blambda = current_blambda / 10
                current_blambda = np.round(current_blambda, precision) # round to 0.1 ** precision
                if current_blambda < lower_blambda or current_blambda in iterated:
                    current_blambda = 10 ** ((np.log10(upper_blambda * lower_blambda)) / 2)
                current_blambda = np.round(current_blambda, precision) # round to 0.1 ** precision
            else:
                lower_blambda = current_blambda
                current_blambda = current_blambda * 10
                current_blambda = np.round(current_blambda, precision) # round to 0.1 ** precision
                if current_blambda > upper_blambda or current_blambda in iterated:
                    current_blambda = 10 ** ((np.log10(upper_blambda * lower_blambda)) / 2)
                current_blambda = np.round(current_blambda, precision) # round to 0.1 ** precision
            iterated.add(current_blambda)
            if (upper_blambda - lower_blambda) / 10 ** (np.floor(np.log10(lower_blambda))) < 3:
                self._output_to_debugfile(f"last cross validate with blambda within interval ({lower_blambda}, {upper_blambda})")
                minor_num = 5
                step = (upper_blambda - lower_blambda) / (minor_num + 2)
                minor_set = np.linspace(lower_blambda + step, upper_blambda - step, minor_num)
                for _ in range(minor_num):
                    current_blambda = minor_set[_]
                    self.blambda_ = current_blambda
                    fitted_nn = self._fit_one_neural_network(i, tr_loader)
                    new_discor = self._dist_cor(fitted_nn, cv_loader)
                    if new_discor > best_dis and self.lower_cov_ <= self.cov_sum_[-1] <= self.upper_cov_:
                        best_dis = new_discor
                        best_blambda = current_blambda
                        self._output_to_debugfile(f"best discor is {best_dis} best blambda is {best_blambda} at last cv cov sum is {self.cov_sum_[-1]}")
                break
        if self.adaptive_cv:
            self._output_to_debugfile(f"blambda adaptive cv has done with best blambda {best_blambda} and best discor {best_dis}")
        current_blambda = best_blambda
        refit_num = 0
        while True:
            ## last tuning to let sample cov in the range while keeping orthogonality.
            self.blambda_ = current_blambda
            if self.early_stop:
                fitted_nn = self._fit_one_neural_network(i, tr_loader, cv_loader)
            else:
                fitted_nn = self._fit_one_neural_network(i, train_loader)

            if self.max_dim == 1 and self.nn_output_dim_ == 1:
                self._output_to_debugfile("no validation for nn with output_dim 1 and max dim 1")
                break

            if self.lower_cov_ <= self.cov_sum_[-1] <= self.upper_cov_:
                if self.method != "seq":
                    pred = self._generate_pred_inbatch(fitted_nn, train_loader)
                    if pred.shape[1] == 1:
                        break
                    normdiff = torch.norm(torch.corrcoef(pred.T) - torch.eye(pred.shape[1], device=self.device)) / pred.shape[1]
                    if normdiff >= 0.316:# sqrt(0.1) ~~ 0.316 
                        # check orthogonality of fnorm penalty
                        self._output_to_debugfile(f"fnorm correlated, {torch.corrcoef(pred.T)[:2, :2]}, diff {normdiff} ,blambda {self.blambda_}")
                        best_blambda = best_blambda * 1.1
                        current_blambda = best_blambda
                    else:
                        self._output_to_debugfile(f"network traininng done with blambda {best_blambda:.2f}, normdiff {normdiff} total refit time is {refit_num} cov sum is {self.cov_sum_[-1]}")
                        break
                elif i >= 1:
                    pred = self._generate_pred_inbatch(fitted_nn, train_loader)
                    for j in range(i):
                        nn = self.estimators_[j]
                        pred = torch.cat([pred, self._generate_pred_inbatch(nn, train_loader)], axis=1)
                    normdiff = torch.norm(torch.corrcoef(pred.T) - torch.eye(pred.shape[1], device=self.device)) / pred.shape[1]
                    self._output_to_debugfile(f"cov {self.cov_sum_[-1]} seq cor nom {normdiff}")
                    if normdiff >= 0.316: # sqrt(0.1) ~~ 0.316 
                        # check orthogonality of fnorm penalty
                        self._output_to_debugfile(f"seq correlated, {torch.corrcoef(pred.T)}, cov {self.cov_sum_[-1]}, tilde {self.blambda_tilde_}, lamb {current_blambda}")
                        self.blambda_tilde_ = self.blambda_tilde_ * 1.1
                    else:
                        self._output_to_debugfile(f"done, {torch.corrcoef(pred.T)}, cov {self.cov_sum_[-1]}, tilde {self.blambda_tilde_}, lamb {current_blambda}")
                        self._output_to_debugfile(f"network traininng done with blambda {best_blambda} total refit time is {refit_num} cov sum is {self.cov_sum_[-1]}")
                        break
                else:
                    self._output_to_debugfile(f"network traininng done, blambda_tilde {self.blambda_tilde_},blambda {best_blambda} total refit time is {refit_num} cov sum is {self.cov_sum_[-1]}")
                    break
            else:
                refit_num += 1
                if refit_num % 2 == 0:
                    self._output_to_debugfile(f"{refit_num} refit, reinital with covsum {self.cov_sum_[-1]}")
                if refit_num % 2 == 0:
                    if self.method == "seq":
                        print(f"cov small {self.cov_sum_[-1]}, tilde {self.blambda_tilde_}, lamb {current_blambda}")
                    else:
                        if self.cov_sum_[-1] < self.lower_cov_:
                            best_blambda = best_blambda / 1.1
                        else:
                            best_blambda = best_blambda * np.log(self.cov_sum_[-1]) / np.log(self.upper_cov_)
                    current_blambda = best_blambda
                if refit_num % 21 == 0:
                    print(f"cov small {self.cov_sum_[-1]}, tilde {self.blambda_tilde_}, lamb {current_blambda}")
                    break
        y_mode_emb_distcor = self._dist_cor(fitted_nn, train_loader)
        if np.isnan(y_mode_emb_distcor) or y_mode_emb_distcor == 0:
            print(y)
            print("y may be the same for all sample")
            # raise Exception("y may be the same for all sample")
        return fitted_nn



    def _fit_one_neural_network(self, i, tr_loader, cv_loader=None):
        ## fit the i th net of estimators
        # if self.early_stop and cv_loader is None:
            # raise Exception("Neural network is fitted with early stopping without specifying cv loader")
        self.gmdd_loss_ = []
        self.total_loss_ = []
        self.cov_loss_ = []
        self.cov_sum_ = []
        self.dc_hist_ = []
        neural_network = copy.deepcopy(self.estimators_[i])
        if self.optimizer == "adam":
            optimizer = torch.optim.Adam(neural_network.parameters(), lr = self.lr)
        elif self.optimizer == "sgd":
            optimizer = torch.optim.SGD(neural_network.parameters(), lr= self.lr, momentum=0.9)
        # optimizer = torch.optim.AdamW(neural_network.parameters(), lr = self.lr)
        best_dc = 0
        best_iternum = 0
        consecutive_decrease_count = 0
        epoch_per_test = 10
        for i_iter in range(self.iter_num):
            self._fit_one_epoch(neural_network, optimizer, i, tr_loader)
            if i_iter % epoch_per_test == 0 and self.early_stop:
                # if i_iter < 100:
                #     continue
                dc = self._dist_cor(neural_network, cv_loader)
                # dc = -self._gmdd_of_loader(neural_network, cv_loader)
                # self._output_to_debugfile(f"{i_iter}th iteration with distance_correlation {dc}")
                self.dc_hist_.append(dc)
                if dc > best_dc:
                    best_dc = dc
                    best_iternum = i_iter
                    consecutive_decrease_count = 0
                    self._output_to_debugfile(f"{i_iter}th iteration with best cv distance_correlation {dc}")
                else:
                    consecutive_decrease_count += 1

            if consecutive_decrease_count >= int(100/epoch_per_test):
                self._output_to_debugfile(f"No distance_correlation improvement on cv loader for consecutive iter at {best_iternum}")
                break

            if self.debug and i_iter % 25 == 0 and self.device != "cpu":
                self._output_to_debugfile(f"{i_iter}th iteration with {torch.cuda.mem_get_info(device=self.device)[0] / 1024**3:2.2f} GB free gpu memory")
        if cv_loader is not None:
            for i_iter in range(best_iternum):
                # self._fit_one_epoch(neural_network, optimizer, i, tr_loader)
                self._fit_one_epoch(neural_network, optimizer, i, cv_loader)
                dc = self._dist_cor(neural_network, cv_loader)
                self.dc_hist_.append(dc)
                # dc = -self._gmdd_of_loader(neural_network, cv_loader)
                # if dc > best_dc:
                    # break
            self._output_to_debugfile(f"early_stop stop at iter {i_iter + 1} with best_iternum {best_iternum}")

        return neural_network

    def _fit_one_epoch(self, neural_network, optimizer, i, tr_loader):
        for temp_x, temp_y in tr_loader:
            yhat = neural_network(temp_x)
            yhatdim = yhat.shape[1]
            tmp_yhat_mean = torch.mean(yhat, axis = 0).reshape(1, -1)
            tmp_centered_yhat = yhat - tmp_yhat_mean.repeat(yhat.shape[0], 1)
            gmdd = self._gmdd_value(yhat, temp_y)
            covsum = torch.sum(torch.abs(torch.cov(tmp_centered_yhat.T)))
            cov_loss = self.blambda_ * self.normfunc_(torch.cov(tmp_centered_yhat.T).reshape(yhatdim, yhatdim) - torch.eye(yhatdim).to(self.device))
            culmulated = 0
            for j in range(i):
                # cov between current net and previous fitted one
                net_ = self.estimators_[j]
                previous = net_(temp_x).detach()
                previous = (previous - torch.mean(previous, axis=0).reshape(1, -1).repeat(temp_x.shape[0], 1)).detach() # centralized
                culmulated += (previous.T.detach() @ tmp_centered_yhat / temp_x.shape[0]) ** 2 # cov  square

            cov_loss = cov_loss + self.blambda_tilde_ * culmulated ## cross cov loss
            loss = cov_loss + gmdd


            self.total_loss_.append(loss.data.cpu().tolist())
            self.gmdd_loss_.append(gmdd.data.cpu().tolist())
            self.cov_loss_.append(cov_loss.data.cpu().tolist())
            self.cov_sum_.append(covsum.data.cpu().tolist())

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

    
    def transform(self, X):
        check_is_fitted(self)
        # X = check_array(X, allow_nd=True)
        X = self._validate_data(X, allow_nd=True, reset=False)
        # construct data loader with second y be empty
        data_loader = torch.utils.data.DataLoader(torch.utils.data.TensorDataset(torch.Tensor(X).to(self.device), ), batch_size=100, shuffle=False) # assure all net get the 

        if self.method == "seq":
            z = torch.zeros((X.shape[0], self.max_dim)).to(self.device)
            for i in range(self.max_dim):
                self.estimators_[i].eval()  # set it to eval 
                z[:, i:(i+1)] = self._generate_pred_inbatch(self.estimators_[i], data_loader)
        else:
            z = self._generate_pred_inbatch(self.estimators_[0], data_loader)
        z = np.array(z.detach().cpu().data, dtype = np.float32)
        
        return z

    def predict(self, X):
        return self.transform(X)
