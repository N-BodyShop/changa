(setq-default c-basic-offset 4)
(defun my-c-mode-common-hook ()
       ;; my customizations for all of c-mode, c++-mode, objc-mode, java-mode
       ;  JST's style
  	(setq c-basic-offset 4)
       (c-set-offset 'block-close '+)
       (c-set-offset 'defun-close '+)
       (c-set-offset 'class-close '+)
       (auto-fill-mode 1)
       (hilit-highlight-buffer)
       )
(add-hook 'c-mode-common-hook 'my-c-mode-common-hook)
